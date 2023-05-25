import argparse
import csv
import requests
import urllib

parser = argparse.ArgumentParser(description="Download AlphaFold PDB files for a list of UniProt accessions.")
parser.add_argument("accession_file", help="Text file containing a list of UniProt accessions.")
parser.add_argument("alphafold_file", help="CSV metadata file containing UniProt accessions with corresponding AlphaFold accessions.")
parser.add_argument("--proteins_directory", default="./proteins", help="Directory to download Uniprot protein FASTA files to. Default is directory named proteins.")
parser.add_argument("--structures_directory", default="./structures", help="Directory to download AlphaFold files to. Default is directory named structures.")

args = parser.parse_args()

# read list of uniprot accessions from file
with open(args.accession_file, "r") as file:
    uniprot_accessions = [line.strip() for line in file]

# Read the AlphaFold file and extract UniProt accessions with their corresponding AlphaFold accessions
accessions_map = {}
with open(args.alphafold_file, "r") as file:
    reader = csv.reader(file)
    for row in reader:
        uniprot_accession = row[0]
        alphafold_accession = row[3]
        accessions_map[uniprot_accession] = alphafold_accession

# Download the AlphaFold files for the UniProt accessions in the input list
for uniprot_accession in uniprot_accessions:
    if uniprot_accession in accessions_map:
        alphafold_accession = accessions_map[uniprot_accession]

        # Download the AlphaFold file
        url = f"https://alphafold.ebi.ac.uk/files/{alphafold_accession}-model_v4.pdb"
        response = requests.get(url)
        pdb_file_path = f"{args.structures_directory}/{uniprot_accession}-F1-model_v1.pdb"
        with open(pdb_file_path, "w") as pdb_file:
            pdb_file.write(response.text)

        print(f"AlphaFold file downloaded for UniProt accession {uniprot_accession} with AlphaFold accession {alphafold_accession}.")

        # Download the UniProt protein FASTA file
        fasta_url = f"https://www.uniprot.org/uniprot/{uniprot_accession}.fasta"
        fasta_file_path = f"{args.proteins_directory}/{uniprot_accession}.fasta"
        urllib.request.urlretrieve(fasta_url, fasta_file_path)
        print(f"UniProt protein FASTA file downloaded for UniProt accession {uniprot_accession}.")

    else:
        print(f"No AlphaFold file found for UniProt accession {uniprot_accession}.")
