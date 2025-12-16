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
import os

def read_tcal_sav(nB, s_type='w', tcal_dir=None, mode='high', date='20190115'):
    """
    nB: int
      beam number
    """
    from scipy.io.idl import readsav
    import os
    if tcal_dir is None:
        tcal_dir = _get_default_tcal_dir()
    fname = os.path.join(tcal_dir, f'{date}/median_{date}.Tcal-results.HI_{s_type}.{mode}.sav')
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
        tcal_dir = _get_default_tcal_dir()
    fname = os.path.join(tcal_dir, f'{date}/CAL.{date}.{mode}.{s_type.upper()}.fits')
    f = fits.open(fname)
    tc_freq = f[1].data['FREQ'][0]
    tc_T = f[1].data['TCAL'][0,nB-1].T
    return tc_freq, tc_T, fname




from .config import conf

def _get_default_tcal_dir():
    raw_dir = conf.get('tcal', 'tcal_dir', '~/Tcal/')
    return os.path.expanduser(raw_dir)

# URLs resolved via config system (Env > Config File > Default)
TCAL_REPO_MANIFEST_URL, TCAL_REPO_BASE_URL = conf.get_tcal_urls()
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
    import os  # Added import os
    from .downloader import download_and_extract_zip

    # 1. Offline Mode Check
    if conf.get_boolean('general', 'offline'):
        return []

    os.makedirs(tcal_dir, exist_ok=True)
    
    timestamp_file = os.path.join(tcal_dir, '.last_update_check')
    manifest_cache_file = os.path.join(tcal_dir, 'manifest_cache.json')
    lock_file_path = os.path.join(tcal_dir, '.update.lock')
    # Phase 1: Separate lock for manifest to avoid race conditions
    manifest_lock_path = os.path.join(tcal_dir, '.manifest.lock') 
    failure_marker = os.path.join(tcal_dir, '.network_failure')
    
    # --- Phase 1: Ensure Manifest Cache is Up-to-Date ---
    # Double-checked locking to prevent redundant downloads
    
    # 1. Optimistic Check
    need_update = True
    if os.path.exists(timestamp_file) and os.path.exists(manifest_cache_file):
        try:
            mtime = os.path.getmtime(timestamp_file)
            if time.time() - mtime < TCAL_UPDATE_INTERVAL:
                need_update = False
        except OSError:
            pass

    if need_update:
        # 2. Acquire Lock and Re-check
        # We use a blocking lock here. Usage of 'a+' matches other lock implementations in this module.
        m_lock = open(manifest_lock_path, 'a+')
        try:
            fcntl.lockf(m_lock, fcntl.LOCK_EX) # Blocking
            
            # Re-check timestamp (someone else might have updated it while we waited)
            already_updated = False
            if os.path.exists(timestamp_file) and os.path.exists(manifest_cache_file):
                try:
                    mtime = os.path.getmtime(timestamp_file)
                    if time.time() - mtime < TCAL_UPDATE_INTERVAL:
                        already_updated = True
                except OSError:
                    pass
            
            if not already_updated:
                try:
                    from .downloader import get_file
                    # Use get_file with overwrite=True to force update
                    res = get_file(TCAL_REPO_MANIFEST_URL, manifest_cache_file, overwrite=True, failure_marker=failure_marker)
                    if res:
                        # Update timestamp
                        with open(timestamp_file, 'w') as f:
                            f.write(str(time.time()))
                except Exception as e:
                    print(f"Warning: Tcal manifest update failed: {e}")
            
        finally:
            fcntl.lockf(m_lock, fcntl.LOCK_UN)
            m_lock.close()

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
            
            # Check for breaker before attempting lock (optimization)
            if os.path.exists(failure_marker):
                 try:
                     if time.time() - os.path.getmtime(failure_marker) < 60:
                          return manifest_dates
                 except: pass

            # If target file is missing, we MUST ensure it gets downloaded.
            # We try to get an exclusive lock.
            lock_file = open(lock_file_path, 'a+') 
            try:
                # First try non-blocking (fast path if we are the only one)
                try:
                    fcntl.lockf(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    # Someone else is holding the lock.
                    # Since we NEED this file and it's missing, we must WAIT.
                    print(f"Waiting for another process to download Tcal data for {target_date}...")
                    fcntl.lockf(lock_file, fcntl.LOCK_EX) # Blocking wait
                
                # Once we have the lock, check existance again (maybe previous owner downloaded it)
                if not os.path.exists(target_path):
                    file_url = f"{TCAL_REPO_BASE_URL}/{target_date}/{target_date}.zip"
                    print(f"Downloading Tcal data for {target_date}...")
                    # Use large timeout (60 retries * 10s = 600s = 10 mins) for slow networks
                    download_and_extract_zip(file_url, target_path, max_retries=60, retry_delay=10, failure_marker=failure_marker)
            
            except Exception as e:
                print(f"Error updating tcal data: {e}")
            finally:
                 try:
                    fcntl.lockf(lock_file, fcntl.LOCK_UN)
                    lock_file.close()
                 except:
                    pass
                 
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
        tcal_dir = _get_default_tcal_dir()
    
    
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

def list_tcal_dates():
    """
    List available Tcal dates with robust fallback:
    1. Online Manifest (Force Update)
    2. Local Cache (if online fails)
    3. Local Directories (if cache fails)
    """
    import sys
    import glob
    import json # Added import json
    tcal_dir = _get_default_tcal_dir()
    
    print(f"Checking Tcal dates in: {tcal_dir}", file=sys.stderr)
    
    dates = []
    source = "online"
    
    # 1. Try Online/Cache via check_and_update_tcal
    # We cheat a bit: passing a non-existent date 'LIST' forces it to load manifest
    # but skip download. However, check_and_update_tcal doesn't force update manifest
    # unless TTL expired or offline.
    # So we manually force update the manifest first.
    
    manifest_url, _ = conf.get_tcal_urls()
    manifest_cache_file = os.path.join(tcal_dir, 'manifest_cache.json')
    failure_marker = os.path.join(tcal_dir, '.network_failure')
    
    try:
        from .downloader import get_file
        if not conf.get_boolean('general', 'offline'):
             res = get_file(manifest_url, manifest_cache_file, overwrite=True, failure_marker=failure_marker)
             if res:
                 source = f"Manifest ({manifest_url})"
             else:
                 source = "Local Cache (Network Failed, using last known good)"
        else:
             print("Warning: Offline mode enabled. List may be incomplete (network unavailable).", file=sys.stderr)
             source = "Local Cache (Offline Mode)"
             
        # Load Cache
        if os.path.exists(manifest_cache_file):
            with open(manifest_cache_file, 'r') as f:
                data = json.load(f)
                dates = data.get('date', [])
                dates.sort()
        else:
            raise FileNotFoundError("Manifest cache not found")

    except Exception as e:
        print(f"Warning: Failed to fetch/read manifest ({e}).", file=sys.stderr)
        source = "local_scan"
        
        # 3. Fallback: Scan Directories
        print("Warning: Falling back to local directory scan.", file=sys.stderr)
        if os.path.exists(tcal_dir):
            # Find directories that look like YYYYMMDD
            dirs = [d for d in os.listdir(tcal_dir) if os.path.isdir(os.path.join(tcal_dir, d)) and d.isdigit() and len(d) == 8]
            dates = sorted(dirs)
    
    if source == "local_scan":
        print(f"Source: Local Scan (Potentially incomplete)", file=sys.stderr)
    else:
        print(f"Source: {source}", file=sys.stderr)
        
    for d in dates:
        print(d)

def download_all_tcal():
    """Downloads all dates available in the manifest."""
    tcal_dir = _get_default_tcal_dir()
    dates = check_and_update_tcal(tcal_dir) # Get manifest
    
    print(f"Found {len(dates)} dates. Starting batch update...")
    for date in dates:
        print(f"Checking/Downloading {date}...")
        check_and_update_tcal(tcal_dir, target_date=date)
    print("All dates processed.")

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description="Manage Tcal Data")
    subparsers = parser.add_subparsers(dest='command')
    
    # list command
    subparsers.add_parser('list', help='List available Tcal dates (Online > Cache > Local)')
    
    # update-all command
    subparsers.add_parser('update-all', help='Download/Update ALL dates from manifest')
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    
    if args.command == 'list':
        list_tcal_dates()
    elif args.command == 'update-all':
        download_all_tcal()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
