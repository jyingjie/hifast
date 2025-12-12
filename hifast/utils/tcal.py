## cal
"""
Tcal Data Management Module
---------------------------

This module handles reading and automatically updating Tcal (temperature calibration) data.

Key Features & Logic:
1. **GitHub-based Storage**: Data is hosted on GitHub Releases. A `manifest.json` lists available dates.
   - Repo: https://github.com/jyingjie/hifast-tcal-data

2. **On-Demand Downloading (Lazy Load)**:
   - Does NOT download all historical data by default.
   - `check_and_update_tcal` fetches the `manifest.json` (cached locally for 24h).
   - Only when a specific date is needed (and missing locally) will it be downloaded.

3. **Smart 'Auto' Selection**:
   - `read_tcal(..., date='auto')` calculates the best match using the UNION of Local files and Remote manifest dates.
   - This ensures the best calibration date is selected even if the file hasn't been downloaded yet.

4. **Concurrency Safety**:
   - Uses `fcntl` non-blocking file locks.
   - Safe to run with `hifast ... -p 50` (multiple parallel processes).
   - Only one process performs the download; others skip or wait.

5. **Offline Mode**:
   - Set environment variable `HIFAST_OFFLINE=1` to disable all network requests.
   - In offline mode, logic degrades gracefully to use ONLY locally available files.
"""
import numpy as np

def read_tcal_sav(nB, s_type='w', tcal_dir=None, mode='high', date='20190115'):
    """
    nB: int
      beam number
    """
    from scipy.io.idl import readsav
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir+f'{date}/median_{date}.Tcal-results.HI_{s_type}.{mode}.sav'
    tc_info= readsav(fname)[f'{mode}_{s_type}'][0]
    tc_freq= tc_info['freq']
    return tc_freq, tc_info['M%02d_TC'%nB], fname

def read_tcal_fits(nB, s_type='w', tcal_dir=None, mode='high', date=''):
    """
    nB: int
      beam number
    """
    from astropy.io import fits
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    fname = tcal_dir + f'{date}/CAL.{date}.{mode}.{s_type.upper()}.fits'
    f = fits.open(fname)
    tc_freq = f[1].data['FREQ'][0]
    tc_T = f[1].data['TCAL'][0,nB-1].T
    return tc_freq, tc_T, fname




# TODO: Replace with actual repository URL
TCAL_REPO_MANIFEST_URL = "https://raw.githubusercontent.com/jyingjie/hifast-tcal-data/main/manifest.json"
TCAL_REPO_BASE_URL = "https://github.com/jyingjie/hifast-tcal-data/releases/download"
TCAL_UPDATE_INTERVAL = 86400  # Check once per 24 hours

def check_and_update_tcal(tcal_dir, target_date=None):
    """
    Checks for Tcal updates, caches the manifest, and optionally downloads a specific date.
    
    Args:
        tcal_dir (str): Local Tcal directory.
        target_date (str, optional): Specific date string (e.g., '20200531') to ensure is available.
        
    Returns:
        list: List of available date strings (remote manifest).
    """
    import time
    import json
    import fcntl
    import requests
    from hifast.utils.downloader import download_and_extract_zip

    # 1. Offline Mode Check
    if os.environ.get('HIFAST_OFFLINE'):
        # In offline mode, try to read cached manifest, otherwise return empty list (or list of local dirs?)
        # Better to return empty list here to avoid blocking, caller should fallback to local glob
        return []

    os.makedirs(tcal_dir, exist_ok=True)
    
    timestamp_file = os.path.join(tcal_dir, '.last_update_check')
    manifest_cache_file = os.path.join(tcal_dir, 'manifest_cache.json')
    lock_file_path = os.path.join(tcal_dir, '.update.lock')
    
    # --- Phase 1: Ensure Manifest Cache is Up-to-Date ---
    need_update = True
    if os.path.exists(timestamp_file):
        try:
            mtime = os.path.getmtime(timestamp_file)
            if time.time() - mtime < TCAL_UPDATE_INTERVAL:
                need_update = False
        except OSError:
            pass
            
    if need_update:
        # Acquire Lock for Manifest Update
        lock_file = open(lock_file_path, 'a+') 
        try:
            fcntl.lockf(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
            
            # Double check inside lock
            if os.path.exists(timestamp_file):
                try:
                    mtime = os.path.getmtime(timestamp_file)
                    if time.time() - mtime < TCAL_UPDATE_INTERVAL:
                        need_update = False
                except OSError:
                    pass
            
            if need_update:
                try:
                    response = requests.get(TCAL_REPO_MANIFEST_URL, timeout=5)
                    if response.status_code == 200:
                        data = response.json()
                        # Save to cache
                        with open(manifest_cache_file, 'w') as f:
                            json.dump(data, f)
                        
                        # Update timestamp
                        with open(timestamp_file, 'w') as f:
                            f.write(str(time.time()))
                except Exception as e:
                    print(f"Warning: Tcal manifest update failed: {e}")
                    
        except BlockingIOError:
            pass # Someone else is updating manifest
        except Exception as e:
            print(f"Error during Tcal update check: {e}")
        finally:
            try:
                fcntl.lockf(lock_file, fcntl.LOCK_UN)
                lock_file.close()
            except:
                pass

    # --- Phase 2: Load Manifest ---
    manifest_dates = []
    if os.path.exists(manifest_cache_file):
        try:
            with open(manifest_cache_file, 'r') as f:
                data = json.load(f)
                manifest_dates = data.get('date', [])
        except:
             pass

    # --- Phase 3: Download Target Date ---
    if target_date and target_date in manifest_dates:
        target_path = os.path.join(tcal_dir, target_date)
        if not os.path.exists(target_path):
            file_url = f"{TCAL_REPO_BASE_URL}/{target_date}/{target_date}.zip"
            # get_file handles its own locking for download
            print(f"Downloading Tcal data for {target_date}...")
            download_and_extract_zip(file_url, target_path)

    return manifest_dates
            
def read_tcal(nB, s_type='w', tcal_dir=None, mode='high', date='auto', mjd=None):
    """
    nB: int
       beam number
    s_type: str
       type, w or n
    tcal_dir: str
       if not set, use ~/Tcal/
    mode: str
       high or low
    date: str
       example: 20190115 or 20200531, default is None
    mjd:
       if date is 'auto', using the nearest date of tcal
    """
    from glob import glob
    import os
    if tcal_dir is None:
        tcal_dir= os.path.expanduser("~")+'/Tcal/'
    
    
    # 1. Update manifest and get list of remote dates
    remote_dates = check_and_update_tcal(tcal_dir) # First pass: just cache manifest

    # 2. Get local dates
    local_dates = [os.path.basename(i) for i in glob(tcal_dir + '/20[0-9][0-9][0-9][0-9][0-9][0-9]')]
    
    # 3. Combine for 'auto' selection
    all_dates = sorted(list(set(remote_dates + local_dates)))

    if len(all_dates) == 0:
         raise(ValueError('can not find tcal file. If this is the first run, ensure internet connection or manually place files in ~/Tcal/'))

    if date == 'auto':
        from astropy.time import Time
        if mjd is None:
            raise(ValueError('need input mjd if date is auto'))
        mjds_h = Time([f'{s[:4]}-{s[4:6]}-{s[6:8]} 00:00:00.000' for s in all_dates], format='iso').mjd
        #Select best date
        date = all_dates[np.argmin(abs(mjds_h - mjd))]
    
    # 4. Trigger download if needed (for both explicit and auto date)
    check_and_update_tcal(tcal_dir, target_date=date)

    # 5. Final check (in case download failed or it's just invalid)
    # We check local glob again ensures we only proceed if file is really there
    # (Optional: optimization to check os.path.exists direct)
    dates_have_now = [os.path.basename(i) for i in glob(tcal_dir + '/20[0-9][0-9][0-9][0-9][0-9][0-9]')]
    
    if date not in dates_have_now:
         raise(ValueError(f'can not find tcal file in {date} (Download may have failed or date invalid)'))
         
    ## read tcal
    if date == '20190115':
        return read_tcal_sav(nB, s_type, tcal_dir, mode, date='20190115')
    else:
        return read_tcal_fits(nB, s_type, tcal_dir, mode, date=date)
