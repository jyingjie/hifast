
import sys
import os
import shutil
import time
import multiprocessing
import traceback

# Add parent directory to path so we can import hifast
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from hifast.utils.tcal import read_tcal, check_and_update_tcal

def worker_process(worker_id, tcal_dir):
    """
    Worker function simulating a single hifast process.
    """
    try:
        # Simulate slight jitter in start times
        time.sleep(worker_id * 0.1) 
        
        # All workers try to read the same file (which is missing initially)
        # This triggers check_and_update_tcal -> locking -> download
        # We use a date we know exists in manifest: 20190115
        print(f"[Worker {worker_id}] Starting read_tcal...")
        
        start_time = time.time()
        # Explicit date to force specific file check/download
        freq, tc, fname = read_tcal(nB=1, tcal_dir=tcal_dir, date='20190115')
        end_time = time.time()
        
        print(f"[Worker {worker_id}] Success! Read {fname} in {end_time - start_time:.2f}s")
        return True
    except Exception as e:
        print(f"[Worker {worker_id}] FAILED: {e}")
        # traceback.print_exc()
        return False

def test_tcal_concurrency():
    print(">>> Testing Tcal Concurrency (Multiprocessing)")
    
    # Use a fresh directory
    test_tcal_dir = os.path.expanduser("~/Tcal_test_concurrency")
    if os.path.exists(test_tcal_dir):
        shutil.rmtree(test_tcal_dir)
    
    print(f"Test Directory: {test_tcal_dir}")
    
    # Reset Environment
    if 'HIFAST_OFFLINE' in os.environ:
        del os.environ['HIFAST_OFFLINE']

    # Number of concurrent processes
    num_workers = 5
    print(f"Spawning {num_workers} processes...")

    pool = multiprocessing.Pool(processes=num_workers)
    
    results = []
    for i in range(num_workers):
        results.append(pool.apply_async(worker_process, (i, test_tcal_dir)))
    
    pool.close()
    pool.join()
    
    # Verify results
    success_count = sum([r.get() for r in results])
    print(f"\n{success_count}/{num_workers} processes succeeded.")
    
    if success_count == num_workers:
        print("PASS: All processes handled concurrency correctly.")
    else:
        print("FAIL: Some processes failed.")

    # Cleanup (Optional, keep if user wants to inspect)
    if os.path.exists(test_tcal_dir):
         shutil.rmtree(test_tcal_dir)
         print(f"Cleanup: Removed {test_tcal_dir}")

if __name__ == "__main__":
    # Ensure multiprocessing works correctly
    multiprocessing.set_start_method('spawn', force=True)
    test_tcal_concurrency()
