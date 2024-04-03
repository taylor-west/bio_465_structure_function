import requests
import re
import os
import pandas as pd
from Bio.KEGG import REST

USING_TEMP_UNIPROT_DATA = True

CURR_DIRECTORY_PATH = os.path.dirname(__file__)
DATAFILES_FILEPATH = os.path.join(os.path.abspath(os.path.join(CURR_DIRECTORY_PATH, os.pardir)), 'datafiles')
TARGET_ORTHOLOGS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_orthologs.csv')
TARGET_ORGANISMS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_organisms.csv')
TARGET_GENES_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_genes.csv')
TEMP_GENES_DIRECTORY= os.path.join(CURR_DIRECTORY_PATH, 'out', 'temp')
TEMP_GENES_FILEPATH = os.path.join(TEMP_GENES_DIRECTORY 'temp_genes.csv')

ORTHOLOGS_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'orthologs_results.csv')
GENES_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'genes_results.csv')
TEMP_UNIPROT_IDS_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'temp', 'temp_uniprot_ids_results.csv')

UNIPROT_IDS_FILEPATH = os.path.join(DATAFILES_FILEPATH, 'orthologs_uniprot_ids.csv')
KO_UNIPROTS_FILEPATH = os.path.join(DATAFILES_FILEPATH, 'ko_uniprots.csv')

# returns a list of uniprot_id values
def find_ortholog_uniprots_by_ko_id(target_organisms_filepath, ko_id):
    print('in kegg_orthologs.find_ortholog_uniprots_by_ko_id...')

    # get genes for the specified KEGG ortholog
    genes = get_genes_by_ortholog(ko_id) # gets gene names for the given ko_id (returns list<string>)
    genes_dict = [{'ko_id': ko_id, 'gene_name': gene_name} for gene_name in genes] # TODO: integrate into line below after testing
    ortholog_genes_df = pd.DataFrame(genes_dict, columns=['ko_id', 'gene_name'])

    # filter the list of genes for the appropriate organisms (if applicable)
    organisms_df = pd.read_csv(target_organisms_filepath) 
    ortholog_genes_df['kegg_organism_code'] = ortholog_genes_df.apply(get_organism_code_from_row, axis=1)
    filtered_genes = ortholog_genes_df.loc[ortholog_genes_df['kegg_organism_code'].isin(organisms_df['kegg_organism_code'])]


    # creates a temporary save file for the uniprots (so you don't have to start completely over if it crashes)
    if not os.path.exists(TEMP_GENES_DIRECTORY):
        os.mkdir(TEMP_GENES_DIRECTORY)
        
    if os.path.exists(TEMP_UNIPROT_IDS_RESULTS_FILEPATH):
        temp_uniprot_df = pd.read_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)
        unprocessed_genes = filtered_genes.loc[filtered_genes['gene_name'].isin(temp_uniprot_df['gene_name']) == False] # only iterate over genes that haven't been already processed
    else:
        temp_uniprot_df = pd.DataFrame(columns=['ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
        unprocessed_genes = filtered_genes

    # get uniprot ids for each unprocessed gene
    for index, row in unprocessed_genes.iterrows():
        uniprot_id = get_uniprot_ids_by_gene(row["gene_name"])

        temp_uniprot_loop_df = pd.DataFrame([{'ko_id': row["ko_id"], 'gene_name': row["gene_name"], 'kegg_organism_code': row["kegg_organism_code"], 'uniprot_id': uniprot_id}], columns=['ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
        temp_uniprot_df = pd.concat([temp_uniprot_df, temp_uniprot_loop_df])
        temp_uniprot_df.to_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH, index=False)


    uniprot_ids_df = temp_uniprot_df.dropna(subset=["uniprot_id"]) 
    organisms_list_filename = os.path.split(target_organisms_filepath)[-1]
    uniprot_ids_df["organisms_list"] = organisms_list_filename  # adds the name of the organisms list file as a new column in the df

    #write uniprot ids to a file
    uniprots_filepath = os.path.join(DATAFILES_FILEPATH, "ortholog_uniprots", f'ko_uniprots_{ko_id}.csv')
    uniprot_ids_df.to_csv(uniprots_filepath)

    if os.path.exists(TEMP_UNIPROT_IDS_RESULTS_FILEPATH):
        os.remove(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)

    uniprot_id_array = uniprot_ids_df["uniprot_id"].array
    return uniprot_id_array

# gets all the orthologs from this KEGG ID
# returns a list of strings
def get_genes_by_ortholog(ko_id):
    print(f'in kegg_orthologs.get_genes_by_ortholog for kegg_ortholog_id={ko_id}...')

    ortholog_genes_response = REST.kegg_find("genes", ko_id).read().rstrip().split('\n')
    
    gene_names = [line.split("\t")[0] for line in ortholog_genes_response]
    return gene_names

def get_orthologs_by_pathway(pathway_name):
    print(f'in kegg_orthologs.get_genes_by_pathway for pathway_name={pathway_name}...')

    url = f'https://rest.kegg.jp/get/{pathway_name}'

    response = requests.get(url)

    if response.status_code == 200:
        pathway_data = response.text
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
    return kegg_ids

def get_uniprot_ids_by_gene(gene_name):
    gene_response = REST.kegg_get(gene_name).read()
    # Regular expression pattern to match UniProt ID
    pattern = r'UniProt:\s+(\S+)'

    # Extracting UniProt ID using regular expression
    match = re.search(pattern, gene_response)
    if match:
        uniprot_id = match.group(1)
    else:
        uniprot_id = None

    # print(f'gene: {gene_name} --> uniprot_id: {uniprot_id}') # TODO: remove after testing
    return uniprot_id

def get_organism_code_from_row(row):
    return row['gene_name'].split(':')[0]

def read_organisms_from_file(filepath):
    pd.read_csv(filepath) 