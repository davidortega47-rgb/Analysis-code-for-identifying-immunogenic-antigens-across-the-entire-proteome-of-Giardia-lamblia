# LIBCE

Laboratorio de Inmunología y Biología Celular – Universidad de Sonora

This repository contains a set of tools designed for identifying immunogenic antigens across the entire proteome of Giardia lamblia.

## Usage 


The iedb_runner.py script automates the MHC-II binding prediction process for protein sequences.

```bash 
python3 iedb_runner.py /path/to/fasta/file
```

## Output Structure

The script generates a LIBCE_RESULTS directory. Each run creates a specific sub-directory following the format: RESULTS-[FastaFileName]-[Species]-[DD-MM-YYYY] (e.g., RESULTS-ProteomaGS-MOUSE-08-12-2025).

Inside the results directory, you will find:

* app.log: A detailed log documenting all processed proteins.

* results_short_[FastaFile].csv: A condensed CSV containing: allele, seq_num, start, end, length, core_peptide, peptide, and ic50.

* results_long_[FastaFile].csv: A comprehensive CSV containing all available data: allele, seq_num, start, end, length, core_peptide, peptide, ic50, rank, adjusted_rank, and protein.

## File Descriptions

* iedb_runner.py: A Python script used to predict MHC-II binding using the NN-align method for Mouse alleles: H2-IAb, H2-IAd, H2-IEd, H2-IAk, and H2-IEk.

* catalog_cleaning.ipynb: A Jupyter notebook for feature engineering of the Proteome Catalog (e.g., ./ProteomaGS/Cat_ProteomaGS.csv). It integrates data to add Cyst, Trophozoite, and Secretome columns.

* file_merge.ipynb: A notebook designed to merge the raw outputs from iedb_runner.py with the cleaned catalog, producing a unified dataset for analysis.

* file_analysis.ipynb: A notebook for Exploratory Data Analysis (EDA) performed on the combined file (e.g., GL50581_Combined.csv) generated in the previous step.

* Requirements.txt: A list of Python libraries required to run this project.

* test_data/: (TBD) Sample data for testing the pipeline.

* GSproteome/: (TBD) Storage for proteome-specific files.

* WBproteome/: (TBD) Storage for proteome-specific files.

## Basic Workflow

The data pipeline follows a sequential flow:

    iedb_runner.py: Generate binding predictions from FASTA files.

    catalog_cleaning.ipynb: Prepare the proteome catalog with specific biological features.

    file_merge.ipynb: Combine prediction results with the biological catalog.

    file_analysis.ipynb: Perform data analysis on the final dataset.

## Installation

Follow these steps to set up the environment and run the LIBCE tools on your local machine.
1. Clone the Repository

Open your terminal and clone the project to your desired directory:

```bash
git clone https://github.com/your-username/LIBCE.git
cd LIBCE
```

2. Set Up a Virtual Environment

It is highly recommended to use a virtual environment to avoid conflicts with other Python packages.

Create the environment:
```bash
    python3 -m venv venv
```

Activate the environment:

        Linux/macOS: source venv/bin/activate

        Windows: .\venv\Scripts\activate

3. Install Dependencies

Once the environment is active, install the required libraries using pip:

```bash
pip install -r Requirements.txt
```


4. External Dependencies (IEDB Tools)

The iedb_runner.py script relies on the IEDB MHC-II binding prediction tool.

Ensure you have the standalone IEDB tools installed locally.
