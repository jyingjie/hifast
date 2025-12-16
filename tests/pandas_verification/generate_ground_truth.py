
import sys
import os
import numpy as np
import pickle
import pandas as pd

# Add hifast-code to path (one level deeper now: tests/pandas_verification/ -> tests/ -> hifast-code/)
# Actually, hifast-code/tests/pandas_verification/.. -> tests/ .. -> hifast-code/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

def generate_gain_para():
    print("--- Generating Gain_para ---")
    try:
        from hifast.core.gain import Gain_para
        df = Gain_para()
        return df
    except Exception as e:
        print(f"Error: {e}")
        return None

def generate_flux_ratio():
    print("--- Generating get_ratio ---")
    try:
        from hifast.core.flux import get_ratio
        ratio1, freq1 = get_ratio(1)
        ratio2, freq2 = get_ratio(1, np.array([1200, 1300]))
        return {'ratio1': ratio1, 'freq1': freq1, 'ratio2': ratio2, 'freq2': freq2}
    except Exception as e:
        print(f"Error: {e}")
        return None

def generate_radec_ky():
    print("--- Generating process_ky ---")
    # Data path relative to this script? tests/data is in ../data/ relative to pandas_verification/
    # verify_pandas_removal used: os.path.join(os.path.dirname(__file__), 'data/KY/...') which was effectively tests/data
    # So now it should be ../data
    ky_file = os.path.join(os.path.dirname(__file__), '../data/KY/M33_OTF_2021_07_31_05_14_00_000.xlsx')
    if not os.path.exists(ky_file):
        print(f"File not found: {ky_file}")
        # Try finding it relative to hifast-code root?
        # tests/data/KY/...
        ky_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../tests/data/KY/M33_OTF_2021_07_31_05_14_00_000.xlsx'))
        
    if not os.path.exists(ky_file):
        print(f"File strictly not found: {ky_file}")
        return None
        
    try:
        from hifast.core.radec import process_ky
        ky_data, mjd, _ = process_ky([ky_file])
        
        clean_data = {}
        for k, v in ky_data.items():
            if hasattr(v, 'values'):
                clean_data[k] = v.values
            else:
                clean_data[k] = v
        
        return {'ky_data': clean_data, 'mjd': mjd}
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    data = {}
    data['gain'] = generate_gain_para()
    data['flux'] = generate_flux_ratio()
    data['radec'] = generate_radec_ky()
    
    out_file = os.path.join(os.path.dirname(__file__), 'ground_truth.pkl')
    with open(out_file, 'wb') as f:
        pickle.dump(data, f)
    print(f"Ground truth saved to {out_file}")
