
import sys
import os
import numpy as np
import re
from glob import glob

# Add hifast-code to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

def test_gain_para():
    print("\n--- Testing Gain_para (hifast.core.gain) ---")
    try:
        from hifast.core.gain import Gain_para
        df = Gain_para()
        print("Type:", type(df))
        # Now returns a dict
        try:
            # data[beam][freq][param]
            val = df['M01'][1050]['a']
            print(f"Sample M01, 'a' at 1050MHz: {val}")
            print(f"Freqs (keys): {list(df['M01'].keys())[:5]} ...")
        except Exception as e:
            print(f"Dict access failed: {e}")

    except Exception as e:
        print("Error testing Gain_para:", e)
        import traceback
        traceback.print_exc()


def test_flux_ratio():
    print("\n--- Testing get_ratio (hifast.core.flux) ---")
    try:
        from hifast.core.flux import get_ratio
        ratio, freq = get_ratio(1)
        print("nB=1 Ratio shape:", ratio.shape)
        print("Freq shape:", freq.shape)
        print("Sample Ratio (first 5):", ratio[:5])
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

def test_radec_process_ky():
    print("\n--- Testing process_ky (hifast.core.radec) ---")
    ky_file = os.path.join(os.path.dirname(__file__), 'data/KY/M33_OTF_2021_07_31_05_14_00_000.xlsx')
    if not os.path.exists(ky_file):
        print(f"Skipping KY test: File not found at {ky_file}")
        return

    from hifast.core.radec import process_ky
    try:
        ky_data, mjd, _ = process_ky([ky_file])
        print(f"KY Data Keys (first 5): {list(ky_data.keys())[:5]}")
        
        # Check specific columns used in radec logic
        cols_to_check = ['SDP_AngleM', 'SDP_PhaPos_X']
        for col in cols_to_check:
            if col in ky_data:
                # Need to handle if it's a pandas Series or a list/array
                val = ky_data[col]
                if hasattr(val, 'values'): # Pandas Series
                     print(f"{col} type: {type(val)} (Pandas)")
                     print(f"{col} sample: {val.values[:3]}")
                else:
                     print(f"{col} type: {type(val)} (Native)")
                     print(f"{col} sample: {np.array(val)[:3]}")
    except Exception as e:
        print("Error processing KY file:", e)
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_gain_para()
    test_flux_ratio()
    test_load_flux_profile()
    test_waterfall_grouping()
    test_radec_process_ky()
