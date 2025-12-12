
import sys
import os
import shutil
import glob
# from astropy.time import Time # Removed to avoid dependency in test runner if possible

# Add parent directory to path so we can import hifast
# File is in tests/utils/, we need to go up two levels to hifast-code/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from hifast.utils.tcal import read_tcal, check_and_update_tcal

def test_tcal_auto_download():
    print(">>> Testing Tcal On-Demand Auto-Update")
    
    # define a test dir to avoid messing with user's real ~/Tcal if possible, 
    # but user said they moved it, so maybe we use default?
    # Let's use a explicit test dir for safety, user can see it works.
    test_tcal_dir = os.path.expanduser("~/Tcal_test_env_rerun")
    if os.path.exists(test_tcal_dir):
        shutil.rmtree(test_tcal_dir)
    
    print(f"Using test directory: {test_tcal_dir}")
    
    # Case 1: Auto date selection (Remote Manifest)
    # Target: 20250329 (Approximately MJD 60763)
    target_mjd = 60763.0 # Hardcoded to avoid using astropy.time.Time here
    print(f"\n[Case 1] Read with date='auto', mjd={target_mjd} (Should match 20250329)")
    
    try:
        # This should:
        # 1. Fetch manifest
        # 2. MATCH 20250329 from remote list
        # 3. DOWNLOAD 20250329
        # 4. Read it
        # Note: read_tcal requires nB (beam number), s_type etc.
        # We assume nB=1 exists in the data.
        freq, tc, fname = read_tcal(nB=1, tcal_dir=test_tcal_dir, date='auto', mjd=target_mjd)
        
        print(f"  Success! Read file: {fname}")
        print(f"  Frequency points: {len(freq)}")
        
        # Verify directory structure
        expected_dir = os.path.join(test_tcal_dir, "20250329")
        if os.path.exists(expected_dir):
            print(f"  Verified: Directory {expected_dir} exists.")
        else:
            print(f"  ERROR: Directory {expected_dir} DOES NOT exist.")
            
    except Exception as e:
        print(f"  FAILED: {e}")
        import traceback
        traceback.print_exc()

    # Case 2: On-Demand caching check
    # Should NOT download fetch manifest again (cached)
    # Should download 20190115 if we ask for it
    print(f"\n[Case 2] Read explicit date='20190115'")
    try:
        freq, tc, fname = read_tcal(nB=1, tcal_dir=test_tcal_dir, date='20190115')
        print(f"  Success! Read file: {fname}")
        
        expected_dir = os.path.join(test_tcal_dir, "20190115")
        if os.path.exists(expected_dir):
             print(f"  Verified: Directory {expected_dir} exists.")
    except Exception as e:
        print(f"  FAILED: {e}")

    # Case 3: Offline Mode Simulation
    print(f"\n[Case 3] Offline Mode (HIFAST_OFFLINE=1)")
    os.environ['HIFAST_OFFLINE'] = '1'
    try:
        # Try to read a date that DOES NOT exist
        # 20200531 is in manifest but not downloaded yet
        # check_and_update_tcal should return empty list or skip download
        # read_tcal should fail saying file not found
        read_tcal(nB=1, tcal_dir=test_tcal_dir, date='20200531')
        print("  ERROR: Should have failed in offline mode for missing file!")
    except ValueError as e:
        print(f"  Expected Failure caught: {e}")
    except Exception as e:
        print(f"  Unexpected error: {e}")
    
    del os.environ['HIFAST_OFFLINE']

    # Case 4: Offline Mode WITH Data
    # We know 20190115 was downloaded in Case 2.
    print(f"\n[Case 4] Offline Mode (HIFAST_OFFLINE=1) WITH Data present")
    os.environ['HIFAST_OFFLINE'] = '1'
    try:
        # A. Explicit read of existing file
        freq, tc, fname = read_tcal(nB=1, tcal_dir=test_tcal_dir, date='20190115')
        print(f"  [A] Success! Read existing file in offline mode: {fname}")
        
        # B. Auto selection locally
        # 20190115 MJD is ~58498. Let's ask for something close to it.
        target_mjd_old = 58500.0 
        freq, tc, fname = read_tcal(nB=1, tcal_dir=test_tcal_dir, date='auto', mjd=target_mjd_old)
        if '20190115' in fname:
             print(f"  [B] Success! Auto-selected local file 20190115 for MJD {target_mjd_old}")
        else:
             print(f"  [B] FAILED: Selected {fname} instead of 20190115")

    except Exception as e:
        print(f"  FAILED: {e}")
        import traceback
        traceback.print_exc()

    del os.environ['HIFAST_OFFLINE']
    
    # Case 5: Online Mode with Non-Existent Directory
    # Should auto-create directory and succeed
    print(f"\n[Case 5] Online Mode with NEW non-existent directory")
    new_tcal_dir = os.path.expanduser("~/Tcal_new_auto_create")
    if os.path.exists(new_tcal_dir):
        shutil.rmtree(new_tcal_dir)
        
    try:
        # Asking for a specific date to trigger download
        # Logic: check_and_update_tcal should create dir, then download manifest, then download file
        freq, tc, fname = read_tcal(nB=1, tcal_dir=new_tcal_dir, date='20190115')
        print(f"  Success! Auto-created directory and downloaded file: {fname}")
        if os.path.exists(new_tcal_dir):
             print(f"  Verified: Directory {new_tcal_dir} was created.")
    except Exception as e:
        print(f"  FAILED: {e}")
        import traceback
        traceback.print_exc()
    finally:
        if os.path.exists(new_tcal_dir):
            shutil.rmtree(new_tcal_dir)

    # Case 6: Offline Mode with Non-Existent Directory
    # Should FAIL with ValueError (Data not found), NOT path error
    print(f"\n[Case 6] Offline Mode with NEW non-existent directory")
    os.environ['HIFAST_OFFLINE'] = '1'
    missing_tcal_dir = os.path.expanduser("~/Tcal_missing_offline")
    if os.path.exists(missing_tcal_dir):
        shutil.rmtree(missing_tcal_dir)
        
    try:
        read_tcal(nB=1, tcal_dir=missing_tcal_dir, date='auto', mjd=50000)
        print("  ERROR: Should have failed!")
    except ValueError as e:
        print(f"  Expected Failure caught: {e}")
        if "can not find tcal file" in str(e):
            print("  Verified: Correct error message.")
    except Exception as e:
        print(f"  Unexpected error: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
    
    del os.environ['HIFAST_OFFLINE']
    print("\nTests Completed.")

if __name__ == "__main__":
    test_tcal_auto_download()
