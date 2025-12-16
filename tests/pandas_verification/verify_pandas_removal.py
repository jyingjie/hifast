
import sys
import os
import numpy as np
import re
from glob import glob
import pickle

# Add hifast-code to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

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
             
             try:
                 beam = 'M01'
                 freq = 1050
                 param = 'a'
                 
                 val_new = data_new[beam][freq][param]
                 
                 # Check if beam column exists in DataFrame
                 if beam in df_old.columns:
                      # df_old[beam] gives a Series where index might satisfy freq?
                      # Or df_old[beam][freq] gives a value (or dict?)
                      # Based on previous debugging: Index was MultiIndex (Param, Freq?) or something.
                      val_old_container = df_old[beam][freq]
                      # val_old_container is the value container.
                      val_old = val_old_container[param]
                      
                     
                      print(f"Old Value: {val_old}, New Value: {val_new}")
                      if np.isclose(float(val_old), float(val_new)):
                          print("SUCCESS: Values match.")
                      else:
                          print(f"FAILURE: Values mismatch: {val_old} vs {val_new}")
                      
                      # Robustness: Check length of keys (freqs)
                      if len(data_new[beam]) != len(df_old.columns): # Assuming df_old cols are just freqs
                          # But wait, df_old cols might be MultiIndex if not parsed? 
                          # df_old is the raw return. 
                          # If it was a DataFrame in pandas version with simple columns (1050, 1100...), check count.
                          # Based on gain.py pandas logic: columns were range(1050, 1500, 50).
                          # Let's verify number of frequency points.
                          n_freqs_new = len(data_new[beam].keys())
                          n_freqs_old = len(df_old.columns)
                          if n_freqs_new == n_freqs_old:
                               print(f"SUCCESS: Frequency count matches ({n_freqs_new}).")
                          else:
                               print(f"FAILURE: Frequency count mismatch: {n_freqs_old} vs {n_freqs_new}")
                 else:
                      # If index is (Beam, Param) and columns are Freq?
                      # Previous error: "Column M01 not found". Index: MultiIndex[(M01, a)...]
                      # So Columns are Frequencies: 1050, 1100...
                      # Access: df.loc[(Beam, Param), Freq]
                      try:
                          val_old = df_old.loc[(beam, param), freq]
                          print(f"Old Value ({beam},{param},{freq}): {val_old}, New: {val_new}")
                          if np.isclose(float(val_old), float(val_new)):
                              print("SUCCESS: Values match.")
                          else:
                              print(f"FAILURE: Values mismatch: {val_old} vs {val_new}")
                      except KeyError:
                          print(f"Key lookup failed in DF structure. Index sample: {df_old.index[:3]} Cols: {df_old.columns[:3]}")
                     
             except Exception as e:
                 print(f"Comparison logic failed: {e}")
        else:
             print("No ground truth for Gain_para.")

        # Basic functional check
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
                if ratio.shape == gt['ratio1'].shape:
                     print(f"SUCCESS: Ratio shape matches {ratio.shape}.")
                else:
                     print(f"FAILURE: Ratio shape mismatch {ratio.shape} vs {gt['ratio1'].shape}.")
            else:
                print("FAILURE: Ratios mismatch.")
                
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
    keys = list(map(lambda x: re.sub('-M[0-1][0-9]','-M00', os.path.basename(x)), fnames))
    
    # Itertools Grouping
    from itertools import groupby
    print("Itertools Grouping Result (Proposed):")
    combined = sorted(zip(keys, fnames), key=lambda x: x[0])
    for key, group in groupby(combined, lambda x: x[0]):
        group_files = [x[1] for x in group]
        print(f"  Key: {key}, Files: {group_files}")

def test_radec_process_ky(ground_truth=None):
    print("\n--- Testing process_ky (hifast.core.radec) ---")
    # path: tests/pandas_verification/ -> tests/data/KY...
    ky_file = os.path.join(os.path.dirname(__file__), '../../tests/data/KY/M33_OTF_2021_07_31_05_14_00_000.xlsx')
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
            keys_new = set(ky_data.keys())
            keys_old = set(gt_data.keys())
            if keys_new == keys_old:
                print("SUCCESS: Keys match exactly.")
            
            cols_to_check = ['SDP_AngleM', 'SDP_PhaPos_X', 'SysTime']
            for col in cols_to_check:
                if col in ky_data and col in gt_data:
                    val_new = np.array(ky_data[col])
                    val_old = np.array(gt_data[col])
                    
                    if np.issubdtype(val_new.dtype, np.number) and np.issubdtype(val_old.dtype, np.number):
                        match = np.allclose(val_new, val_old, equal_nan=True)
                    else:
                        match = np.array_equal(val_new, val_old)
                        
                    if match:
                        print(f"SUCCESS: Column {col} matches.")
                    else:
                        print(f"FAILURE: Column {col} mismatch.")
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
    test_load_flux_profile()
    test_waterfall_grouping()
    test_radec_process_ky(gt)
