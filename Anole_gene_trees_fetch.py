import requests
import csv
import json
import os
import time
import signal
from tqdm import tqdm
import timeout_decorator
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Global variables to track progress
last_processed_gene = None
processed_genes = set()
total_genes = 0
current_gene_number = 0

# Function to fetch gene information from Ensembl with timeout
@timeout_decorator.timeout(300)  # 5 minutes timeout
def fetch_gene_info(gene_id):
    lookup_url = f"https://rest.ensembl.org/lookup/id/{gene_id}?content-type=application/json"
    response = requests.get(lookup_url, verify=False)
    print(f"Fetching gene info: Status Code {response.status_code}")
    print(f"URL used: {lookup_url}")
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch gene information. Response: {response.text}")
        return None

# Function to fetch gene tree information from Ensembl with timeout
@timeout_decorator.timeout(300)  # 5 minutes timeout
def fetch_gene_tree_info(gene_id, gene_symbol):
    gene_tree_url = f"https://rest.ensembl.org/genetree/member/symbol/anolis_carolinensis/{gene_symbol}?content-type=application/json"
    print(f"Fetching gene tree info. URL: {gene_tree_url}")
    response = requests.get(gene_tree_url, verify=False)
    print(f"Fetching gene tree info: Status Code {response.status_code}")
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch gene tree information. Response: {response.text}")
        return None

# Function to process gene tree data
def process_gene_tree_data(gene_tree_info):
    processed_data = []

    def traverse_tree(node):
        if 'children' in node:
            for child in node['children']:
                traverse_tree(child)
        else:
            if 'taxonomy' in node:
                species_name = node['taxonomy'].get('scientific_name', 'N/A')
                gene_id = node.get('id', 'N/A')
                gene_name = node.get('gene_member', {}).get('display_name', 'N/A')
                processed_data.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'species': species_name,
                })

    if isinstance(gene_tree_info, list) and len(gene_tree_info) > 0:
        if 'tree' in gene_tree_info[0]:
            traverse_tree(gene_tree_info[0]['tree'])
    elif isinstance(gene_tree_info, dict) and 'tree' in gene_tree_info:
        traverse_tree(gene_tree_info['tree'])

    return processed_data

# Function to save checkpoint
def save_checkpoint(processed_genes):
    with open('checkpoint.json', 'w') as f:
        json.dump({
            'processed_genes': list(processed_genes),
            'last_gene': last_processed_gene,
            'current_gene_number': current_gene_number
        }, f)
    print(f"Checkpoint saved. Last processed gene: {last_processed_gene}")

# Function to load checkpoint
def load_checkpoint():
    if os.path.exists('checkpoint.json'):
        with open('checkpoint.json', 'r') as f:
            checkpoint_data = json.load(f)
            return (
                set(checkpoint_data['processed_genes']),
                checkpoint_data['last_gene'],
                checkpoint_data.get('current_gene_number', 0)
            )
    return set(), None, 0

# Signal handler function
def signal_handler(signum, frame):
    print("\nProcess interrupted. Saving checkpoint...")
    save_checkpoint(processed_genes)
    print("You can resume later by running the script again.")
    exit(0)

# Function to process a batch of genes
def process_gene_batch(batch, output_dir, total_genes):
    global last_processed_gene, processed_genes, current_gene_number

    for gene in tqdm(batch, desc="Processing genes", unit="gene"):
        current_gene_number += 1
        # Skip if gene has already been processed
        if gene['gene_id'] in processed_genes:
            print(f"Skipping already processed gene: {gene['gene_id']}")
            continue

        print(f"\nProcessing gene {current_gene_number} of {total_genes}: {gene['gene_id']}")

        try:
            # Fetch gene information
            gene_info = fetch_gene_info(gene['gene_id'])
            if not gene_info:
                print(f"Failed to retrieve gene information for {gene['gene_id']}. Skipping.")
                continue

            gene_symbol = gene_info.get('display_name', gene['gene_symbol'])

            # Fetch gene tree information
            gene_tree_info = fetch_gene_tree_info(gene['gene_id'], gene_symbol)

            # Process gene tree data or write "No gene tree available"
            if gene_tree_info:
                processed_data = process_gene_tree_data(gene_tree_info)

                if processed_data:
                    # Write processed data to CSV
                    output_file = os.path.join(output_dir, f'{gene_symbol}_gene_tree.csv')
                    with open(output_file, 'w', newline='') as csvfile:
                        fieldnames = ['gene_id', 'gene_name', 'species']
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                        for row in processed_data:
                            writer.writerow(row)

                    print(f"Gene tree information for {gene_symbol} has been written to {output_file}")
                    print(f"Number of entries: {len(processed_data)}")
                else:
                    print(f"No gene tree data found for {gene_symbol} after processing.")
                    output_file = os.path.join(output_dir, f'{gene_symbol}_gene_tree.txt')
                    with open(output_file, 'w') as txtfile:
                        txtfile.write("No gene tree available")
                    print(f"No gene tree available for {gene_symbol}. Written to {output_file}")
            else:
                output_file = os.path.join(output_dir, f'{gene_symbol}_gene_tree.txt')
                with open(output_file, 'w') as txtfile:
                    txtfile.write("No gene tree available")
                print(f"No gene tree available for {gene_symbol}. Written to {output_file}")

            print(f"Successfully processed {gene_symbol}")

            # Add processed gene to the set and update last processed gene
            processed_genes.add(gene['gene_id'])
            last_processed_gene = gene['gene_id']

            # Save checkpoint after each gene
            save_checkpoint(processed_genes)

        except timeout_decorator.TimeoutError:
            print(f"Timeout occurred while processing gene {gene['gene_id']}. Moving to next gene.")
        except Exception as e:
            print(f"An error occurred while processing gene {gene['gene_id']}: {str(e)}")

        # Optional: Add a small delay to prevent overwhelming the API
        time.sleep(1)

# Main function to process gene tree information for all Green Anole protein-coding genes
def process_all_gene_trees():
    global last_processed_gene, processed_genes, total_genes, current_gene_number

    # Create a directory to store all gene tree files
    output_dir = 'green_anole_gene_tree_files'
    os.makedirs(output_dir, exist_ok=True)

    # Load checkpoint if exists
    processed_genes, last_gene, current_gene_number = load_checkpoint()

    # Read the Green Anole protein-coding genes CSV file
    anole_genes = []
    with open('Green_Anole_protein-coding_genes.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
        print(f"CSV columns: {fieldnames}")

        for i, row in enumerate(reader):
            if i < 5:  # Print first 5 rows for debugging
                print(f"Row {i+1}: {row}")

            gene_id = row.get('gene_id') or row.get('ensembl_id') or row.get('WBGeneID')
            gene_symbol = row.get('gene_symbol') or row.get('symbol') or gene_id

            if not gene_id:
                print(f"Warning: No gene ID found for row: {row}")
                continue

            anole_genes.append({
                'gene_id': gene_id,
                'gene_symbol': gene_symbol
            })

    total_genes = len(anole_genes)
    print(f"Total Green Anole protein-coding genes: {total_genes}")
    print(f"Starting from gene number: {current_gene_number + 1}")

    # Process genes in batches of 100
    batch_size = 100
    for i in range(current_gene_number, len(anole_genes), batch_size):
        batch = anole_genes[i:i+batch_size]
        print(f"\nProcessing batch {i//batch_size + 1} of {len(anole_genes)//batch_size + 1}")
        process_gene_batch(batch, output_dir, total_genes)

    print("\nAll genes have been processed.")
    save_checkpoint(processed_genes)

# Run the main function
if __name__ == "__main__":
    # Set up signal handler
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        process_all_gene_trees()
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Saving checkpoint before exiting...")
        save_checkpoint(processed_genes)
        print("You can resume later by running the script again.")
