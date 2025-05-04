import requests
import csv

def fetch_green_anole_genes():
    # Using the US East Ensembl mirror
    biomart_url = "https://useast.ensembl.org/biomart/martservice"

    xml_query = """<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
        <Dataset name="acarolinensis_gene_ensembl" interface="default">
            <Attribute name="ensembl_gene_id"/>
            <Attribute name="external_gene_name"/>
            <Attribute name="gene_biotype"/>
        </Dataset>
    </Query>
    """

    try:
        response = requests.get(biomart_url, params={'query': xml_query})
        response.raise_for_status()

        lines = response.text.strip().split('\n')

        if len(lines) < 2:
            print(f"Unexpected response format. Response content:\n{response.text}")
            return None

        header = lines[0].split('\t')
        print(f"Response header: {header}")

        genes = []
        for line in lines[1:]:
            parts = line.split('\t')
            if len(parts) != 3:
                print(f"Unexpected line format: {line}")
                continue
            gene_id, gene_name, biotype = parts
            if biotype == 'protein_coding':
                genes.append({
                    'ensembl_id': gene_id,
                    'gene_symbol': gene_name,
                    'biotype': biotype
                })

        print(f"Fetched {len(genes)} protein-coding genes.")
        return genes

    except requests.exceptions.RequestException as e:
        print(f"Error fetching genes: {e}")
        if 'response' in locals():
            print(f"Response content: {response.text}")
        return None

def create_gene_csv(genes):
    with open('Green_Anole_protein-coding_genes.csv', 'w', newline='') as csvfile:
        fieldnames = ['gene_id', 'gene_symbol', 'ensembl_id']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for gene in genes:
            writer.writerow({
                'gene_id': gene['ensembl_id'],  # Using Ensembl ID as gene_id
                'gene_symbol': gene['gene_symbol'],
                'ensembl_id': gene['ensembl_id']
            })

    print(f"CSV file 'Green_Anole_protein-coding_genes.csv' has been created with {len(genes)} genes.")

# Fetch genes and create CSV
genes = fetch_green_anole_genes()
if genes:
    create_gene_csv(genes)
else:
    print("Failed to create CSV file due to error in fetching genes.")
