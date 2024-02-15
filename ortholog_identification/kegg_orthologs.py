import requests
import re
import os
# import json
import pandas as pd
from Bio.KEGG import REST

PATHWAY_KEYWORDS = ["metabolic", "metabolism"]
CURR_DIRECTORY_PATH = os.path.dirname(__file__)
READ_PATHWAY_FROM_FILE = True
READ_ORTHOLOGS_FROM_FILE = True
READ_GENES_FROM_FILE = False

# Step 1: get a list of the correct pathways for humans
# Step 2: find all genes/entries associated with those pathways
# Step 3: get KEGG ortholog entry from gene
# Step 4: get all genes associated with that ortholog, making sure to track organisms

def get_uniprot_ids_for_pathways():
    print('in kegg_orthologs.get_uniprot_ids_for_pathways...')

    # get pathways data
    if READ_PATHWAY_FROM_FILE:
       pathways_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'target_pathways.csv')) 
       
    else:
        pathways = get_pathways()
        pathways_df = pd.DataFrame(pathways, columns=['pathway_name', 'pathway_description'])

        # write pathways to a file (FOR TESTING/DEV)
        pathways_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'target_pathways.csv'))

    # get orthologs for each pathway
    if READ_ORTHOLOGS_FROM_FILE:
       pathway_orthologs_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'pathway_orthologs.csv')) 
       
    else:
        kegg_ids_list = []
        for index, row in pathways_df.iterrows():
            kegg_ids = get_orthologs_by_pathway(row['pathway_name'])
            kegg_ids_list.append(kegg_ids)
        pathway_orthologs_df = pathways_df.assign(ko_id=kegg_ids_list).explode('ko_id').dropna(subset=["ko_id"])       

        # write orthologs to a file (FOR TESTING/DEV)
        pathway_orthologs_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'pathway_orthologs.csv'))




    # get genes for each KEGG ortholog
    if READ_GENES_FROM_FILE:
        ortholog_genes_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'ortholog_genes.csv')) 
    else:
        # ortholog_genes_df = df = pathway_orthologs_df.reindex(columns=["pathway_name", 'pathway_description', "ko_id", "gene_name", "gene_organism"])
        gene_list = []
        for index, row in pathway_orthologs_df.iterrows():
            genes = get_genes_by_ortholog(row["ko_id"]) # gets an array of the gene names
            gene_list.append(genes)
            # ortholog_genes_df.merge(genes, on='ko_id', how='left') # old way using a merge, but is super slow
        ortholog_genes_df = pathway_orthologs_df.assign(gene_name=gene_list).explode('gene_name').dropna(subset=["gene_name"])       
        
        #write genes to a file
        ortholog_genes_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'ortholog_genes.csv'))


    # get uniprot ids for each KEGG ortholog
    uniprot_ids = []
    for index, row in ortholog_genes_df.iloc[:, "gene_name"]:
        uniprot_ids.append(get_uniprot_ids_by_gene(row["gene_name"]))
    uniprot_ids_df = ortholog_genes_df.assign(uniprot_id=uniprot_ids)
    

    #write uniprot ids to a file
    uniprot_ids_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'uniprot_ids.csv'))
    
    ortholog_data = []
    return ortholog_data

# retrieves a list of metabolic pathways matching the PATHWAY_KEYWORDS
def get_pathways():
    print('in kegg_orthologs.get_pathways...')

    # Gets a list of all human pathways
    human_pathways = REST.kegg_list("pathway", "hsa").read().rstrip().split('\n')

    target_pathways = []
    print(f"{len(human_pathways)} pathways found.")

    # Filter human pathways for target keywords
    for line in human_pathways:
        entry, description = line.split("\t")
        if any(keyword.lower() in description.lower() for keyword in PATHWAY_KEYWORDS):
            target_pathways.append([entry, description])

    print(f"{len(target_pathways)} pathways matched the keywords.")

    # converts the pathway names and descriptions into a Pandas Dataframe
    target_pathways_df = pd.DataFrame(target_pathways, columns=['pathway_name', 'pathway_description'])

    return target_pathways_df


# gets all the orthologs from this KEGG ID
def get_genes_by_ortholog(kegg_ortholog_id):
    print(f'in kegg_orthologs.get_genes_by_ortholog for kegg_ortholog_id={kegg_ortholog_id}...')

    ortholog_genes_response = REST.kegg_find("genes", kegg_ortholog_id).read().rstrip().split('\n')
    
    ortholog_genes_response_df = pd.DataFrame([line.split("\t") for line in ortholog_genes_response], columns=['gene_name', 'gene_description'])

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
    print('in kegg_orthologs.get_uniprot_ids_by_gene...')
    gene_response = REST.kegg_get(gene_name).read()
    # TODO: extract uniprot_id from KEGG gene response
    gene_response
    return

if __name__ == "__main__":
    get_uniprot_ids_for_pathways()


# Exactly what we want to do, but it is in R instead of Python
# https://support.bioconductor.org/p/109871/