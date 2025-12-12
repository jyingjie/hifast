

__all__ = ['download_to_temp_and_move', 'get_file', 'download_and_extract_zip']


import os
import shutil
import tempfile
import requests
import fcntl  # For file locking
import time   # For sleep during retry

def download_to_temp_and_move(url, save_path, resume=False):
    """
    Downloads a file from a URL to a temporary file, then moves it to the final location
    upon successful completion.

    Args:
        url (str): The URL of the file to download.
        save_path (str): The full path where the file should be saved.
        resume (bool, optional): If True, try to resume an existing partial download. Defaults to False.

    Returns:
        bool: True if the download was successful (or resumed successfully), False otherwise.
    """

    temp_dir = os.path.dirname(save_path)
    with tempfile.NamedTemporaryFile(delete=False, prefix=os.path.basename(save_path)+'.', dir=temp_dir) as temp_file:
        temp_path = temp_file.name

        headers = {}
        existing_size = 0
        if resume and os.path.exists(temp_path):
            existing_size = os.path.getsize(temp_path)
            headers['Range'] = f'bytes={existing_size}-'

        try:
            response = requests.get(url, stream=True, headers=headers)

            if response.status_code == 200:
                total_size = int(response.headers.get('content-length', 0))
            elif response.status_code == 206 and resume:
                content_range = response.headers.get('content-range')
                if content_range:
                    total_size = int(content_range.split('/')[-1])
                else:
                    total_size = 0
            else:
                print(f"Download failed: status code {response.status_code}")
                return None

            os.makedirs(os.path.dirname(temp_path), exist_ok=True)

            with open(temp_path, 'ab' if resume else 'wb') as file:
                downloaded_size = existing_size
                for chunk in response.iter_content(chunk_size=1048576):
                    file.write(chunk)
                    downloaded_size += len(chunk)

                    # Progress indication (optional)
                    print(f"\rDownloaded {downloaded_size}/{total_size} bytes...", end='')

            print(f"\nFile downloaded successfully to: {temp_path}")  # Print temp path

            # Move the temporary file to the final destination
            shutil.move(temp_path, save_path)
            print(f"File moved to: {save_path}")

            return save_path

        except requests.exceptions.RequestException as e:
            print(f"Download error: {e}")
            print(f"Failed to download file from: {url}, you may need to manually download it and put it at {save_path}.")
            if not resume:
                try:
                    os.remove(temp_path)
                except OSError:
                    pass
            return None


def get_file(url, save_path, resume=False, max_retries=6, retry_delay=10, overwrite=False, failure_marker=None):
    # --- Circuit Breaker Check ---
    if failure_marker and os.path.exists(failure_marker):
        try:
            # If failure was less than 60 seconds ago, skip
            if time.time() - os.path.getmtime(failure_marker) < 60:
                # print(f"Skipping download due to recent network failure: {failure_marker}")
                return None
        except OSError:
            pass

    # check if file already exists
    if not overwrite and os.path.exists(save_path):
        print(f"File already exists at: {save_path}")
        return save_path
    # mkdir if dirname doesn't exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    lock_file_path = save_path + ".lock"
    itsme = True
    for attempt in range(max_retries):
        lock_file = open(lock_file_path, 'w')
        try:
            fcntl.lockf(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)

            if not overwrite and os.path.exists(save_path):
                print(f"File already exists at: {save_path}")
                return save_path
            else:
                if itsme or overwrite:
                    res = download_to_temp_and_move(url, save_path, resume)
                    
                    # --- Circuit Breaker Update ---
                    if failure_marker:
                        try:
                            if res:
                                if os.path.exists(failure_marker):
                                    os.remove(failure_marker)
                            else:
                                with open(failure_marker, 'w') as f:
                                    f.write(str(time.time()))
                        except OSError:
                            pass
                    return res
                else:
                    # Should not reach here if non-blocking lock acquired
                    print(f"Another process has tried to download '{save_path}' and failed.",
                          f"You may need to manually download it from {url} and put it at {save_path}.")
                    return None
        except BlockingIOError:  # Lock is held by another process
            itsme = False
            if attempt < max_retries - 1:
                # If we are waiting for a lock, it means someone else is downloading.
                # If we wanted to OVERWRITE, we might still want to do it after they finish? 
                # Or if we just want "latest", getting their version is fine?
                # For manifest: if someone else is updating, their update is fresh enough. We can accept it.
                # So if locking conflict occurs, we treat it as "someone else handled it".
                # BUT, `get_file` retries. When it gets the lock eventually, `os.path.exists` will be true.
                # If `overwrite=True`, we would download AGAIN.
                # This ensures we get OUR version? No, for manifest, if A updates, B waiting -> B gets lock -> B updates again.
                # This is redundant but consistent with "overwrite=True".
                if attempt == 0:
                    print(f"Another process is downloading '{save_path}'. Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                print(f"Another process is downloading '{save_path} and max retries reached.'. Skipping...")
                return None
        finally:
            try:
                # check if the file can be locked by me
                # only remove lock file if it's created by me
                fcntl.lockf(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.lockf(lock_file, fcntl.LOCK_UN) # Release the lock
                lock_file.close()
                os.remove(lock_file_path)
            except:
                pass

def download_and_extract_zip(url, target_dir, resume=False, max_retries=6, retry_delay=10, failure_marker=None):
    """
    Downloads a ZIP file and extracts it to the target directory atomically.
    Ensures that the target directory is only created if extraction is successful.
    
    Args:
        url (str): The URL of the ZIP file.
        target_dir (str): The directory where the contents should be extracted.
        resume (bool): Whether to resume download.
        max_retries (int): Number of retries for lock acquisition.
        retry_delay (int): Delay in seconds between retries.
        failure_marker (str, optional): Path to circuit breaker file.
        
    Returns:
        str or None: path to target_dir if successful, None otherwise.
    """
    if os.path.exists(target_dir):
        return target_dir

    # 1. Download ZIP to a temporary location
    # Use get_file to handle locking and downloading safely
    # We append .zip to the directory name to create a temp zip path
    # But get_file expects a target path.
    # We'll use a temp directory for the zip file
    
    zip_target_path = os.path.join(os.path.dirname(target_dir), f".tmp_{os.path.basename(target_dir)}.zip")
    
    downloaded_zip = get_file(url, zip_target_path, resume=resume, max_retries=max_retries, retry_delay=retry_delay, failure_marker=failure_marker)
    
    if not downloaded_zip:
        return None
        
    # 2. Extract to a temp dir
    temp_extract_dir = os.path.join(os.path.dirname(target_dir), f".tmp_extract_{os.path.basename(target_dir)}")
    if os.path.exists(temp_extract_dir):
         shutil.rmtree(temp_extract_dir)
         
    try:
        import zipfile
        with zipfile.ZipFile(downloaded_zip, 'r') as zip_ref:
            zip_ref.extractall(temp_extract_dir)
            
        # 3. Atomic Move
        # If target_dir appeared in the meantime (race condition), we should check again
        if os.path.exists(target_dir):
             shutil.rmtree(temp_extract_dir)
             return target_dir
             
        os.rename(temp_extract_dir, target_dir)
        
        # Cleanup zip file
        try:
            os.remove(downloaded_zip)
        except OSError:
            pass
            
        return target_dir
        
    except Exception as e:
        print(f"Error extracting zip: {e}")
        if os.path.exists(temp_extract_dir):
            shutil.rmtree(temp_extract_dir)
        return None
