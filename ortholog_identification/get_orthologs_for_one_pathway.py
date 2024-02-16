import requests
import re
import os
# import json
import pandas as pd
from Bio.KEGG import REST
import sys

PATHWAY_KEYWORDS = ["metabolic", "metabolism"]
CURR_DIRECTORY_PATH = os.path.dirname(__file__)
READ_GENES_FROM_FILE = False


# Step 1: get a list of the correct pathways for humans
# Step 2: find all genes/entries associated with those pathways
# Step 3: get KEGG ortholog entry from gene
# Step 4: get all genes associated with that ortholog, making sure to track organisms

def get_uniprot_ids_for_pathway(pathway_id):
    print('in kegg_orthologs.get_uniprot_ids_for_pathway...')

    ortholog_df = pd.DataFrame(columns=['pathway_id', 'kegg_id', 'gene_id', 'uniprot_id'])

    # get KEGG id for the pathway
    kegg_ids = get_orthologs_by_pathway(pathway_id)

    # get genes for the KEGG ortholog id
    for kegg_id in kegg_ids:
        genes = get_genes_by_ortholog(kegg_id)  # gets an array of the gene names

        # get uniprot ids for each gene
        uniprot_ids = []
        for gene in genes:
            uniprot_id = get_uniprot_ids_by_gene(gene)
            uniprot_ids.append(uniprot_id)
            ortholog_df.loc[len(ortholog_df)] = [pathway_id, kegg_id, gene, uniprot_id]

    print(ortholog_df)
    return ortholog_df


# gets all the orthologs from this KEGG ID
def get_genes_by_ortholog(kegg_ortholog_id):
    print(f'in kegg_orthologs.get_genes_by_ortholog for kegg_ortholog_id={kegg_ortholog_id}...')

    ortholog_genes_response = REST.kegg_find("genes", kegg_ortholog_id).read().rstrip().split('\n')

    ortholog_genes_response_df = pd.DataFrame([line.split("\t") for line in ortholog_genes_response],
                                              columns=['gene_name', 'gene_description'])

    return ortholog_genes_response_df['gene_name']


def get_orthologs_by_pathway(pathway_name):
    # https://www.biostars.org/p/6224/#6355
    print(f'in kegg_orthologs.get_genes_by_pathway for pathway_name={pathway_name}...')

    url = f'https://rest.kegg.jp/get/{pathway_name}'

    response = requests.get(url)

    if response.status_code == 200:
        pathway_data = response.text
        # print(pathway_data)
    else:
        return None

    kegg_ids = []

    # Split the pathway data into lines
    lines = pathway_data.split('\n')
    # Process each line
    in_gene_section = False
    for line in lines:
        # Process each line as needed
        if line.startswith('GENE'):
            in_gene_section = True

        if in_gene_section == False:
            continue
        else:
            if line.startswith('COMPOUND'):
                break
            # extract the KO id using regex
            # Define the regular expression pattern to match 'KO:K12047'
            pattern = r'\[KO:(K\d+)\]'
            # Search for the pattern in the line
            match = re.search(pattern, line)
            # If a match is found, extract the desired information
            if match:
                ko_code = match.group(1)
                kegg_ids.append(ko_code)
                print("Extracted KO code:", ko_code)
            else:
                print("KO code not found in the line.")

    return kegg_ids
    # see below link for some information on Bio.KEGG.REST
    # https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html
    # kegg_get()


def get_uniprot_ids_by_gene(gene_name):
    print(f'in kegg_orthologs.get_uniprot_ids_by_gene {gene_name}...')
    gene_response = REST.kegg_get(gene_name).read()
    # Regular expression pattern to match UniProt ID
    pattern = r'UniProt:\s+(\S+)'
    # Extracting UniProt ID using regular expression
    match = re.search(pattern, gene_response)
    if match:
        uniprot_id = match.group(1)
        # print("UniProt ID:", uniprot_id)
    else:
        uniprot_id = None
        # print("UniProt ID not found.")
    return uniprot_id


if __name__ == "__main__":
    args = sys.argv[1:]
    pathway_id = args[0]
    get_uniprot_ids_for_pathway(pathway_id)

# Exactly what we want to do, but it is in R instead of Python
# https://support.bioconductor.org/p/109871/