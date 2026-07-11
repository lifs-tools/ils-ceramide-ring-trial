#!/usr/bin/env python3
"""
Script to generate 8 separate JSON datasets from the CSV file using the Summarized Study Result Schema.
Each dataset contains data for one matrix (SRM, hTAG, DB, YAA) and one data type (All, Filt).

Usage:
    python3 generate_datasets.py

Requirements:
    - CSV file at: manuscript/output/Suppl_TableS01_Concentrations-CVs_ALL-vs-Outlierfilt-datasets.csv
    - Schema: summarized-study-result-schema.json

Output:
    - 8 JSON files in the datasets/ directory
"""

import csv
import json
from pathlib import Path

# CSV file path
CSV_FILE = Path("manuscript/output/Suppl_TableS01_Concentrations-CVs_ALL-vs-Outlierfilt-datasets.csv")
OUTPUT_DIR = Path("datasets")
OUTPUT_DIR.mkdir(exist_ok=True)

# Publication metadata from CrossRef
PUBLICATION_METADATA = {
    "title": "Concordant inter-laboratory derived concentrations of ceramides in human plasma reference materials via authentic standards",
    "authors": [
        "Federico Torta", "Nils Hoffmann", "Bo Burla", "Irina Alecu", "Makoto Arita",
        "Takeshi Bamba", "Steffany A. L. Bennett", "Justine Bertrand-Michel", "Britta Br\u00fcgger",
        "M\u00f3nica P. Cala", "Dolores Camacho-Mu\u00f1oz", "Antonio Checa", "Michael Chen",
        "Michaela Chocholou\u0161kov\u00e1", "Michelle Cinel", "Emeline Chu-Van", "Benoit Colsch",
        "Cristina Coman", "Lisa Connell", "Bebiana C. Sousa", "Alex M. Dickens",
        "Maria Fedorova", "Finnur Freyr Eir\u00edksson", "Hector Gallart-Ayala", "Mohan Ghorasaini",
        "Martin Giera", "Xue Li Guan", "Mark Haid", "Thomas Hankemeier", "Amy Harms",
        "Marcus H\u00f6ring", "Michal Hol\u010dapek", "Thorsten Hornemann", "Chunxiu Hu",
        "Andreas J. H\u00fclsmeier", "Kevin Huynh", "Christina M. Jones", "Julijana Ivanisevic",
        "Yoshihiro Izumi", "Harald C. K\u00f6feler", "Sin Man Lam", "Mike Lange",
        "Jong Cheol Lee", "Gerhard Liebisch", "Katrice Lippa", "Andrea F. Lopez-Clavijo",
        "Malena Manzi", "Manuela R. Martinefski", "Raviswamy G. H. Math", "Satyajit Mayor",
        "Peter J. Meikle", "Mar\u00eda Eugenia Monge", "Myeong Hee Moon", "Sneha Muralidharan",
        "Anna Nicolaou", "Thao Nguyen-Tran", "Valerie B. O'Donnell", "Matej Ore\u0161i\u010d",
        "Arvind Ramanathan", "Fabien Riols", "Daisuke Saigusa", "Tracey B. Schock",
        "Heidi Schwartz-Zimmermann", "Guanghou Shui", "Madhulika Singh",
        "Masatomo Takahashi", "Margr\u00e9t Thorsteinsd\u00f3ttir", "Noriyuki Tomiyasu",
        "Anthony Tournadre", "Hiroshi Tsugawa", "Victoria J. Tyrrell",
        "Grace van der Gugten", "Michael O. Wakelam", "Craig E. Wheelock",
        "Denise Wolrab", "Guowang Xu", "Tianrun Xu", "John A. Bowden",
        "Kim Ekroos", "Robert Ahrends", "Markus R. Wenk"
    ],
    "year": 2024,
    "doi": "10.1038/s41467-024-52087-x",
    "description": "Nature Communications publication of the interlaboratory ceramide ring trial study"
}

STUDY_DESIGN = {
    "accession": "OBI:0003658",
    "name": "study design type",
    "value": "ring trial design",
    "cvLabel": "OBI"
}

QUANTITY_UNIT = {
    "accession": "UO:0000276",
    "name": "micromole per liter",
    "cvLabel": "UO"
}

# Controlled vocabularies used in the datasets
CONTROLLED_VOCABULARIES = [
    {
        "name": "Ontology for Biomedical Investigations",
        "uri": "https://obi-ontology.org/",
        "version": "2024-07-11"
    },
    {
        "name": "Unit Ontology",
        "uri": "https://www.ebi.ac.uk/ontology-units/",
        "version": "2024-07-11"
    }
]

# Material abbreviations mapping
MATERIAL_NAMES = {
    "DB": "NIST RM 8231 Diabetic (DB)",
    "YAA": "NIST RM 8231 Young African-American (YAA)",
    "T2D": "NIST RM 8231 Type-2 Diabetic (T2D)",
    "SRM": "NIST SRM 1950 (SRM)",
    "hTAG": "hTAG"
}


def read_csv_data():
    """Read CSV file and return list of data rows."""
    rows = []
    
    with open(CSV_FILE, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        
        # Skip header rows
        header_complete = False
        header_parts = []
        
        for row in reader:
            if not row or all(cell.strip() == '' for cell in row):
                continue
            
            # Check if this is still a header row
            if not header_complete:
                # Collect all header parts
                header_parts.extend([cell.strip().strip('"') for cell in row if cell.strip()])
                # Check if we have all expected columns (we expect at least Matrix, Ceramide, etc.)
                if any('Matrix' in cell for cell in row):
                    continue
                # If we hit a row that doesn't contain header keywords, it's a data row
                if not any(kw in ''.join(row) for kw in ['Matrix', 'Ceramide', 'n (', '\u00b5mol/L', '%CV', '%RCV']):
                    # This is actually a data row, so add the header we collected
                    header_complete = True
                    # Put the row back to be processed as data
                    # We need to process the header parts we collected
                    pass
            
            if header_complete:
                # Clean and combine the row
                cleaned_row = []
                for cell in row:
                    cell = cell.strip().strip('"')
                    if cell:
                        cleaned_row.append(cell)
                if cleaned_row:
                    rows.append(cleaned_row)
    
    return rows


def create_dataset(matrix, data_type, rows):
    """Create a dataset JSON structure using the Summarized Study Result Schema."""
    
    lipid_quantities = []
    
    # Get full matrix name
    full_matrix_name = MATERIAL_NAMES.get(matrix, matrix)
    
    for row in rows:
        if len(row) < 14:
            continue
        
        row_matrix = row[0]
        if row_matrix != matrix:
            continue
        
        lipid = row[1]
        
        # Determine which columns to use based on data_type
        if data_type == "All":
            n = int(float(row[2]))
            mean = float(row[3])
            sd = float(row[4])
            median = float(row[5])
            cv = float(row[6])
            rcv = float(row[7])
            col_offset = 2
        else:  # Filt
            n = int(float(row[8]))
            mean = float(row[9])
            sd = float(row[10])
            median = float(row[11])
            cv = float(row[12])
            rcv = float(row[13])
            col_offset = 8
        
        # Determine dataType value
        data_type_value = "Without Outliers" if data_type == "Filt" else "All"
        
        quantity = {
            "lipids": [lipid],
            "quantityUnit": QUANTITY_UNIT,
            "groupingAttributes": {
                "sampleMatrix": {
                    "accession": "OBI:0000747",
                    "name": "material sample",
                    "value": full_matrix_name,
                    "cvLabel": "OBI"
                },
                "dataType": {
                    "accession": "OBI:0001755",
                    "name": "selection criterion",
                    "value": data_type_value,
                    "cvLabel": "OBI"
                }
            },
            "stats": {
                "n": n,
                "mean": mean,
                "sd": sd,
                "median": median,
                "cv": cv,
                "rcv": rcv
            }
        }
        
        lipid_quantities.append(quantity)
    
    # Build description
    if data_type == "All":
        desc = f"Consensus concentration values across all labs for {full_matrix_name} matrix samples from the interlaboratory comparison study"
    else:
        desc = f"Consensus concentration values with outliers removed for {full_matrix_name} matrix samples from the interlaboratory comparison study. Outliers were removed based on Tukey's 1.5 \u00d7 IQR fences, represented and determined separately for each of the four ceramides and reference materials."
    
    dataset = {
        "title": f"Ceramide Concentrations - {full_matrix_name} Matrix - {data_type} Labs",
        "description": desc,
        "nativeId": f"ILS_Ceramide_RingTrial_{matrix}_{data_type}",
        "controlledVocabularies": CONTROLLED_VOCABULARIES,
        "studyDesign": STUDY_DESIGN,
        "publicationMetadata": PUBLICATION_METADATA,
        "lipidSummarizedQuantities": lipid_quantities,
        "visibility": "PUBLIC"
    }
    
    return dataset


def main():
    # Read CSV
    rows = read_csv_data()
    print(f"Read {len(rows)} data rows")
    
    # Get unique matrices
    matrices = sorted(set(row[0] for row in rows if row))
    print(f"Matrices: {matrices}")
    
    # Generate 8 datasets (4 matrices \u00d7 2 types)
    for matrix in matrices:
        for data_type in ["All", "Filt"]:
            dataset = create_dataset(matrix, data_type, rows)
            
            # Generate filename
            matrix_key = matrix.lower().replace(" ", "_")
            filename = OUTPUT_DIR / f"ring_trial_{matrix_key}_{data_type.lower()}.json"
            
            with open(filename, 'w') as f:
                json.dump(dataset, f, indent=2)
            
            print(f"Created: {filename}")
            print(f"  - {len(dataset['lipidSummarizedQuantities'])} lipid entries")


if __name__ == "__main__":
    main()
