#!/usr/bin/env python3
"""
Скрипт проверяет наличие обновленной версии current_tree.json в репозитории YFullTeam.
Использует HTTP ETag (If-None-Match) для предотвращения скачивания без необходимости.
"""
import os
import urllib.request
import urllib.error
import logging

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

URL = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TREE_FILE = os.path.join(SCRIPT_DIR, "yleaf", "data", "hg_prediction_tables", "tree.json")
ETAG_FILE = os.path.join(SCRIPT_DIR, ".tree.etag")

def check_and_update():
    logging.info("Checking for YFull tree updates...")
    
    local_etag = ""
    if os.path.exists(ETAG_FILE):
        with open(ETAG_FILE, "r") as f:
            local_etag = f.read().strip()
            
    req = urllib.request.Request(URL)
    if local_etag:
        req.add_header("If-None-Match", local_etag)
        
    try:
        with urllib.request.urlopen(req) as response:
            remote_etag = response.getheader("ETag")
            
            logging.info("New tree version found! Downloading ~3MB...")
            data = response.read()
            
            # Validate JSON minimally before overwriting
            if not data.strip().startswith(b"{"):
                logging.error("Downloaded data does not look like valid JSON! Aborting.")
                return

            with open(TREE_FILE, "wb") as f:
                f.write(data)
                
            if remote_etag:
                with open(ETAG_FILE, "w") as f:
                    f.write(remote_etag)
                    
            logging.info(f"Successfully updated tree.json (ETag: {remote_etag})")
            
    except urllib.error.HTTPError as e:
        if e.code == 304:
            logging.info("Local YFull tree is already up to date. No download needed.")
        else:
            logging.error(f"HTTP Error checking for updates: {e.code} - {e.reason}")
    except Exception as e:
        logging.error(f"Failed to check/update YTree: {e}")

if __name__ == "__main__":
    check_and_update()
