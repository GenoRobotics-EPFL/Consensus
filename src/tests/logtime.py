import logging
import time
import os
import datetime

class LogTime:
    def __init__(self, description: str):
        self.description = description

    def __enter__(self):
        self.start_time = time.time()

    def __exit__(self, exc_type, exc_value, traceback):
        elapsed_time = time.time() - self.start_time
        logging.info(f'{self.description} took {elapsed_time:.2f} seconds.')

def initialize_log_file(filename: str):
    log_dir = os.path.join("assets", "logs")
    base_filename = os.path.basename(filename).split('.')[0]  # Extract name without extension
    log_filepath = os.path.join(log_dir, f"{base_filename}_computation_time.log")

    # Create a unique filename if the log file already exists
    counter = 1
    while os.path.exists(log_filepath):
        log_filepath = os.path.join(log_dir, f"{base_filename}_computation_time_{counter}.log")
        counter += 1

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logging.basicConfig(filename=log_filepath, level=logging.INFO, format='%(asctime)s - %(message)s')

    # File details
    file_size_bytes = os.path.getsize(filename)
    file_size_mb = file_size_bytes / 1_048_576  # convert bytes to megabytes
    with open(filename, 'r') as f:
        num_lines = sum(1 for line in f)
    
    # File timestamps
    creation_time = datetime.datetime.fromtimestamp(os.path.getctime(filename))
    modified_time = datetime.datetime.fromtimestamp(os.path.getmtime(filename))

    # Add file information at the beginning
    with open(log_filepath, 'a') as log_file:
        log_file.write(f"Processing file: {filename}\n")
        log_file.write(f"File size: {file_size_mb:.2f} MB\n")
        log_file.write(f"Number of lines: {num_lines}\n")
        log_file.write(f"Created on: {creation_time}\n")
        log_file.write(f"Last modified on: {modified_time}\n")
        log_file.write(f"Full path: {os.path.abspath(filename)}\n")
        log_file.write("=" * 50 + "\n")