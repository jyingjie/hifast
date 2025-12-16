# Pandas Removal Verification

This directory contains scripts to verify that replacing `pandas` with standard libraries in `hifast` produces identical results.

## Files

-   `verify_pandas_removal.py`: The main test script. It runs specific checks on `Gain_para`, `get_ratio`, and `process_ky`. If `ground_truth.pkl` exists, it performs a strict comparison against that baseline.
-   `generate_ground_truth.py`: A script to generate the `ground_truth.pkl` file. **This must be run in an environment/commit where `pandas` is still used.**
-   `ground_truth.pkl`: The baseline data (not committed to git, generate it locally).

## How to Run Verification

### 1. Generate Ground Truth (Baseline)

You need to switch to the commit where `pandas` was still part of the codebase (e.g., the state of `v1.4` before removal).

**Commit Hash**: `98bcab2d606164ee422b62172e7d34183d3046d1`

```bash
# 1. Stash any current changes
git stash

# 2. Checkout the pandas-dependent version
git checkout 98bcab2d606164ee422b62172e7d34183d3046d1

# 3. Generate the ground truth file
python tests/pandas_verification/generate_ground_truth.py
# Reference output: Ground truth saved to .../tests/pandas_verification/ground_truth.pkl
```

### 2. Verify Refactored Code

Now switch back to your feature branch (e.g., `feature/remove-pandas` or `v1.4` after merge).

```bash
# 1. Switch back
git checkout feature/remove-pandas
# or
git switch -

# 2. Run the verification
python tests/pandas_verification/verify_pandas_removal.py
```

Expected Output:
-   `Gain_para`: SUCCESS: Values match.
-   `get_ratio`: SUCCESS: Ratios match.
-   `process_ky`: SUCCESS: Keys match exactly... Column matches.

### Notes
-   The `ground_truth.pkl` file is excluded from git to avoid storing large binary blobs.
-   `generate_ground_truth.py` is kept to allow regeneration if logic changes in the future (though you'd need the old code logic to compare against).
