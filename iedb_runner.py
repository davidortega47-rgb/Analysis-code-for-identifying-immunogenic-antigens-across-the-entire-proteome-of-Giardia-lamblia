import sys
import time
import logging
from datetime import date
from pathlib import Path
from multiprocessing import Pool
from random import randint
from time import sleep
from typing import List, Dict, Union
from Bio import SeqIO
import iedb
import pandas as pd


# Config 
NUM_PROCESSES = 4
SPECIES_TARGET = 'mouse'
MAX_API_RETRIES = 15     
RETRY_DELAY_SECONDS = 10
METHOD = "nn_align-2.3"

ALLELES: Dict[str, str] = {
    'human': """DPA1*01/DPB1*04:01,DPA1*01:03/DPB1*02:01,DPA1*02:01/DPB1*01:01,
                DPA1*02:01/DPB1*05:01,DPA1*03:01/DPB1*04:02,DQA1*01:01/DQB1*05:01,
                DQA1*01:02/DQB1*06:02,DQA1*03:01/DQB1*03:02,DQA1*04:01/DQB1*04:02,
                DQA1*05:01/DQB1*02:01,DQA1*05:01/DQB1*03:01,DRB1*01:01,DRB1*03:01,
                DRB1*04:01,DRB1*04:05,DRB1*07:01,DRB1*08:02,DRB1*09:01,DRB1*11:01,
                DRB1*12:01,DRB1*13:02,DRB1*15:01,DRB3*01:01,DRB3*02:02,DRB4*01:01,
                DRB5*01:01""", 
    'mouse': """H2-IAb,H2-IAd,H2-IEd,H2-IAk,H2-IEk"""
}

MOUSE_ALLELES_STR = ALLELES[SPECIES_TARGET]


# Path handling
try:
    FASTA_PATH = Path(sys.argv[1])
except IndexError:
    print("Error: Please provide the path to the FASTA file as a command-line argument.")
    sys.exit(1)

# Path component extraction
FASTA_FILENAME = FASTA_PATH.name
FASTA_NAME = FASTA_FILENAME.split('.')[0]
CURRENT_DATE = date.today().strftime("%d-%m-%Y")

# Define all paths
ROOT_RESULTS = Path(f"../LIBCE_RESULTS/RESULTS-{FASTA_NAME}-{SPECIES_TARGET.upper()}-{CURRENT_DATE}")
LONG_DIR_PATH = ROOT_RESULTS / f"results_long_{FASTA_NAME}"
SHORT_DIR_PATH = ROOT_RESULTS / f"results_short_{FASTA_NAME}"
LOG_PATH = ROOT_RESULTS / "app.log"

# Create all directories
try:
    for directory in [ROOT_RESULTS, LONG_DIR_PATH, SHORT_DIR_PATH]:
        directory.mkdir(parents=True, exist_ok=True)
    print(f"Results directories created under: {ROOT_RESULTS.resolve()}")
except OSError as e:
    print(f"Error creating directories: {e}")
    sys.exit(1)

# Configure Logging
logging.basicConfig(
    filename=LOG_PATH, filemode='w',
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S'
)
logging.info(f'Starting process for FASTA file: {FASTA_PATH.resolve()}')

# Retry mechanism
def query_mhcii_with_retry(sequence: str, alleles_str: str, max_retries: int, delay_seconds: int, protein_name: str) -> Union[pd.DataFrame, None]:
    """
    Handles the IEDB API call with a retry mechanism for transient failures.
    Retries if the call raises an exception OR if the call returns None.
    Returns DataFrame on success, None on final failure.
    """
    attempts = 0
    while attempts < max_retries:
        attempts += 1
        
        sleep(randint(1, 6)) 
        
        try:
            # API call
            mhcii_res = iedb.query_mhcii_binding(
                method=METHOD, 
                sequence=sequence, 
                allele=alleles_str, 
                length=None
            )
            
  
            if mhcii_res is not None and not mhcii_res.empty:
                logging.info(f'{protein_name} - IEDB query successful on attempt {attempts}.')
                return mhcii_res
            elif mhcii_res is not None and mhcii_res.empty:
                 logging.warning(f'{protein_name} - IEDB query successful, but returned an EMPTY DataFrame. Not retrying, treating as completion.')
                 return mhcii_res 
            else:
                 raise Exception("API call returned None unexpectedly.")
            
        # Catch exceptions (e.g., network errors, timeouts, status codes from iedb)
        except Exception as e:
            if attempts == max_retries:
                logging.error(f'{protein_name} - Max retries reached ({max_retries}). Failing this sequence. Last Error: {e}')
                return None 

            logging.warning(f'{protein_name} - IEDB query failed (Attempt {attempts}/{max_retries}): {e}')
            logging.info(f'{protein_name} - Retrying in {delay_seconds} seconds...')
            time.sleep(delay_seconds)
    
    return None



def iedb_worker(seq_record, long_dir: Path, short_dir: Path, alleles_str: str, species: str):
    """
    Worker function to process a single sequence record, using the robust retry function.
    """
    sleep(randint(1, 10)) 
    
    protein_name, sequence = seq_record.id, str(seq_record.seq).replace('X', 'A').replace('B', 'A').replace('J', 'A')
    print(f'Processing: {protein_name}')
    
    # Run IEDB binding prediction using the retry function
    mhcii_res = query_mhcii_with_retry(
        sequence=sequence,
        alleles_str=alleles_str,
        max_retries=MAX_API_RETRIES,
        delay_seconds=RETRY_DELAY_SECONDS,
        protein_name=protein_name
    )
    
    if mhcii_res is None:
         return 

    
    

    # Save results
    # 1. Generating long DataFrame (full results)
    logging.info(f'{protein_name} Writing files.')
    long_result_path = long_dir / f"{protein_name}_RESULTS_LONG.csv"
    mhcii_res['protein'] = protein_name
    mhcii_res.to_csv(long_result_path, index=False)
    
    # 2. Generating short DataFrame (truncated results)
    short_result_path = short_dir / f"{protein_name}_RESULTS_SHORT.csv"
    mhcii_res.iloc[:, :8].to_csv(short_result_path, index=False)
    
    logging.info(f'{protein_name} Has been completed and saved.')




def main():
    """Main execution function for parallel processing."""
    
    try:
        fasta_sequences = list(SeqIO.parse(FASTA_PATH, 'fasta'))
        seq_number = len(fasta_sequences)
        logging.info(f'Number of sequences: {seq_number}')
    except Exception as e:
        logging.error(f"Error reading FASTA file: {e}")
        sys.exit(1)

    if seq_number == 0:
        logging.warning("No sequences found in the FASTA file. Exiting.")
        return

    config_tuple = (LONG_DIR_PATH, SHORT_DIR_PATH, MOUSE_ALLELES_STR, SPECIES_TARGET)
    
    tasks_iterable = [(seq, *config_tuple) for seq in fasta_sequences]
    
    start_time = time.time()
    logging.info('Process Start')

    try:
        with Pool(NUM_PROCESSES) as p:
            p.starmap(iedb_worker, tasks_iterable)
            
    except Exception as e:
        logging.critical(f"A CRITICAL error occurred during multiprocessing Pool setup: {e}", exc_info=True)

    end_time = time.time()
    total_time = end_time - start_time
    
    logging.info('Process End')
    logging.info(f'Total time taken: {total_time:.2f} seconds.')
    print(f'Done! Results saved under: {ROOT_RESULTS.resolve()}')


if __name__ == '__main__':
    main()