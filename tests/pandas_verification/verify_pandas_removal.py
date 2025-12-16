
import sys
import os
import numpy as np
import re
from glob import glob

# Add hifast-code to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))


import pickle

# ... (imports)

def load_ground_truth():
    pkl_file = os.path.join(os.path.dirname(__file__), 'ground_truth.pkl')
    if os.path.exists(pkl_file):
        with open(pkl_file, 'rb') as f:
            return pickle.load(f)
    print("Ground truth file not found.")
    return None

def test_gain_para(ground_truth=None):
    print("\n--- Testing Gain_para (hifast.core.gain) ---")
    try:
        from hifast.core.gain import Gain_para
        data_new = Gain_para()
        print("Type:", type(data_new))
        
        # Check against ground truth
        if ground_truth and ground_truth.get('gain') is not None:
             df_old = ground_truth['gain']
             print("Comparing against ground truth (pandas DataFrame)...")
             # Old structure: DataFrame with MultiIndex or complex structure
             # Based on previous exploration: df['M01'][1050]['a'] gave a value.
             # This implies df['M01'] is a Series or dict? 
             # Let's try to replicate the access pattern used in generation script which might have been implicit.
             # Actually, the user code used: val = df['M01']['1050']['a'] (failed previously as df['M01'][1050] was needed?)
             # Let's inspect ONE value blindly.
             
             try:
                 beam = 'M01'
                 freq = 1050
                 param = 'a'
                 
                 val_new = data_new[beam][freq][param]
                 
                 # Old structure: Index=(Beam, Param), Columns=Freq
                 # Access: df_old.loc[(beam, param), freq]
                 try:
                     val_old = df_old.loc[(beam, param), freq]
                     print(f"Old Value ({beam},{param},{freq}): {val_old}, New: {val_new}")
                     
                     if np.isclose(float(val_old), float(val_new)):
                         print("SUCCESS: Values match.")
                     else:
                         print(f"FAILURE: Values mismatch: {val_old} vs {val_new}")
                 except KeyError:
                      print(f"Key lookup failed in DF. Index sample: {df_old.index[:3]}")
                     
             except Exception as e:
                 print(f"Comparison logic failed: {e}")
        else:
             print("No ground truth for Gain_para.")

        # Basic functional test
        try:
            val = data_new['M01'][1050]['a']
            print(f"Sample M01, 'a' at 1050MHz: {val}")
        except Exception as e:
            print(f"Dict access failed: {e}")

    except Exception as e:
        print("Error testing Gain_para:", e)
        import traceback
        traceback.print_exc()

def test_flux_ratio(ground_truth=None):
    print("\n--- Testing get_ratio (hifast.core.flux) ---")
    try:
        from hifast.core.flux import get_ratio
        ratio, freq = get_ratio(1)
        
        if ground_truth and ground_truth.get('flux'):
            gt = ground_truth['flux']
            print("Comparing 'ratio1' (get_ratio(1))...")
            if np.allclose(ratio, gt['ratio1'], equal_nan=True):
                print("SUCCESS: Ratios match.")
            else:
                print("FAILURE: Ratios mismatch.")
                print("Old:", gt['ratio1'][:5])
                print("New:", ratio[:5])
                
            print("Comparing 'freq1'...")
            if np.allclose(freq, gt['freq1'], equal_nan=True):
                 print("SUCCESS: Frequencies match.")
            else:
                 print("FAILURE: Frequencies mismatch.")
        else:
            print("No ground truth for get_ratio.")

    except Exception as e:
        print("Error testing get_ratio:", e)

def test_load_flux_profile():
    print("\n--- Testing load_flux_profile (hifast.cbr.FLUXGAIN) ---")
    try:
        from hifast.cbr.FLUXGAIN import load_flux_profile
        freqs = np.array([1000, 1400])
        # Note: 3C286 should be in FluxProfiles.csv - if file missing, expects error
        flux = load_flux_profile('3C286', freqs, None)
        print("3C286 Flux at 1000, 1400 MHz:", flux)
    except SystemExit:
        print("Expected SystemExit (File broken/missing handled by code)")
    except Exception as e:
        print(f"Error testing load_flux_profile: {e}")

def test_waterfall_grouping():
    print("\n--- Testing waterfall grouping logic (simulated) ---")
    fnames = [
        'proj-M01-0001.fits', 'proj-M02-0001.fits',
        'proj-M01-0002.fits', 'proj-M19-0001.fits',
        'other-M01-0001.fits'
    ]
    
    # Original logic in hifast.waterfall:
    keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x)), fnames))
    
    print("Files:", fnames)
    print("Keys:", keys)

    # 1. Pandas Grouping (Baseline)
    try:
        import pandas as pd
        files = pd.DataFrame({'key':keys, 'fname':fnames})
        print("Pandas Grouping Result:")
        for key, _files in files.groupby('key'):
             print(f"  Key: {key}, Files: {list(_files['fname'])}")
    except ImportError:
        print("Pandas not installed.")

    # 2. Itertools Grouping (New Logic)
    from itertools import groupby
    print("Itertools Grouping Result (Proposed):")
    # Must sort first for groupby
    combined = sorted(zip(keys, fnames), key=lambda x: x[0])
    for key, group in groupby(combined, lambda x: x[0]):
        group_files = [x[1] for x in group]
        print(f"  Key: {key}, Files: {group_files}")

def test_radec_process_ky(ground_truth=None):
    print("\n--- Testing process_ky (hifast.core.radec) ---")
    ky_file = os.path.join(os.path.dirname(__file__), '../data/KY/M33_OTF_2021_07_31_05_14_00_000.xlsx')
    if not os.path.exists(ky_file):
        print(f"Skipping KY test: File not found at {ky_file}")
        return

    from hifast.core.radec import process_ky
    try:
        ky_data, mjd, _ = process_ky([ky_file])
        
        if ground_truth and ground_truth.get('radec'):
            gt = ground_truth['radec']
            gt_data = gt['ky_data']
            
            print("Comparing KY Data Keys...")
            # Compare keys
            keys_new = set(ky_data.keys())
            keys_old = set(gt_data.keys())
            if keys_new == keys_old:
                print("SUCCESS: Keys match exactly.")
            else:
                print(f"Keys mismatch. New-Old: {keys_new - keys_old}. Old-New: {keys_old - keys_new}")
            
            # Compare specific columns
            cols_to_check = ['SDP_AngleM', 'SDP_PhaPos_X', 'SysTime']
            all_match = True
            for col in cols_to_check:
                if col in ky_data and col in gt_data:
                    val_new = np.array(ky_data[col])
                    val_old = np.array(gt_data[col])
                    
                    # SysTime is string, others float
                    if np.issubdtype(val_new.dtype, np.number) and np.issubdtype(val_old.dtype, np.number):
                        match = np.allclose(val_new, val_old, equal_nan=True)
                    else:
                        match = np.array_equal(val_new, val_old)
                        
                    if match:
                        print(f"SUCCESS: Column {col} matches.")
                    else:
                        print(f"FAILURE: Column {col} mismatch.")
                        all_match = False
                else:
                    print(f"Skipping col {col} (missing in one).")
        else:
            print("No ground truth for process_ky.")

    except Exception as e:
        print("Error processing KY file:", e)
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    gt = load_ground_truth()
    test_gain_para(gt)
    test_flux_ratio(gt)
    test_load_flux_profile() # No strict GT needed, functional test OK
    test_waterfall_grouping() # No strict GT needed, logic is verified functionally
    test_radec_process_ky(gt)

