import os
import scanpy as sc
import datetime
import pandas as pd
import argparse
from dotenv import load_dotenv
load_dotenv()

print("=== Script started ===")

# Fix Python path for importing from scripts/
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), "scripts"))

from cleaning import clean_adata, annotate_adata
from ontology import map_cell_ontology
from metadata import add_metadata

# Argument parser for cleaning and annotation methods
parser = argparse.ArgumentParser(description="Clean h5ad files with different methods.")
parser.add_argument("--cleaning", type=str, default="code", choices=["code", "openai", "claude"], help="Cleaning method to use.")
parser.add_argument("--annotation", type=str, default="code", choices=["code", "openai"], help="Annotation method to use.")
args = parser.parse_args()

print("=== About to process files in data/ ===")

# Paths
data_folder = "data"
output_folder = os.path.join(data_folder, "cleaned")
os.makedirs(output_folder, exist_ok=True)

log_path = os.path.join(output_folder, "cleaning_log.txt")
log_entries = []

# Process all .h5ad files in data/
for filename in os.listdir(data_folder):
    print(f"Found file: {filename}")
    if filename.endswith(".h5ad") and not filename.startswith("cleaned_"):
        filepath = os.path.join(data_folder, filename)
        print(f"🧼 Cleaning {filename} with cleaning='{args.cleaning}', annotation='{args.annotation}' ...")

        # Load
        adata = sc.read_h5ad(filepath)
        before_shape = adata.shape

        # Clean
        adata = clean_adata(adata, method=args.cleaning)
        # Annotate
        adata = annotate_adata(adata, method=args.annotation, name=filename)
        adata = map_cell_ontology(adata)
        adata = add_metadata(adata, name=filename, source_url="https://cellxgene.cziscience.com", description="Cleaned dataset from hackathon")

        # Save
        if args.cleaning == "openai":
            out_file = os.path.join(output_folder, f"openai_{filename}")
        elif args.cleaning == "claude":
            out_file = os.path.join(output_folder, f"claude_{filename}")
        else:
            out_file = os.path.join(output_folder, f"cleaned_{filename}")
        adata.write(out_file)
        after_shape = adata.shape

        # Log
        log_entry = (
            f"{filename}\n"
            f" - Cells: {before_shape[0]} → {after_shape[0]}\n"
            f" - Genes: {before_shape[1]} → {after_shape[1]}\n"
            f" - Ontologies: {adata.obs.get('cell_type_ontology', pd.Series()).nunique() if 'cell_type_ontology' in adata.obs else 'N/A'}\n"
            f"--------------------------------------\n"
        )
        log_entries.append(log_entry)

# Write the full log
explanation = (
    "# Cleaning Log Explanation\n"
    "Each block below corresponds to a single .h5ad file that was cleaned.\n"
    "For each file, you will see:\n"
    "- The filename of the original dataset.\n"
    "- The number of cells (rows) before and after cleaning.\n"
    "- The number of genes (columns) before and after cleaning.\n"
    "- The number of unique cell type ontology annotations present in the cleaned data.\n"
    "\nExample:\n"
    "c4ff5829-a4ef-411e-8469-c884e494eaa0.h5ad\n"
    " - Cells: 17528 → 17528\n"
    " - Genes: 33137 → 19507\n"
    " - Ontologies: 1\n"
    "--------------------------------------\n"
    "\n"
    "- If the number of cells or genes drops, it means the cleaning process removed some based on quality or expression criteria.\n"
    "- The number of unique ontologies gives you a sense of how many distinct cell types were annotated in the cleaned data.\n"
    "\n"
)
with open(log_path, "w") as f:
    f.write(explanation)
    f.writelines(log_entries)

print("✅ All datasets cleaned and saved to 'data/cleaned/'")
print("📄 Log written to:", log_path)
