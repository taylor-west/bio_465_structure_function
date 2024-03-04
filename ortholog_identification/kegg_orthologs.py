import requests
import re
import os
import pandas as pd
from Bio.KEGG import REST

# READ_PATHWAYS_FROM_FILE = True
# READ_ORTHOLOGS_FROM_FILE = True
# READ_GENES_FROM_FILE = True
# READ_ORGANISMS_FROM_FILE = True
# USING_TEMP_GENE_DATA = True
USING_TEMP_UNIPROT_DATA = True
# PATHWAY_KEYWORDS = ["metabolic", "metabolism"]

CURR_DIRECTORY_PATH = os.path.dirname(__file__)
DATAFILES_FILEPATH = os.path.join(os.path.abspath(os.path.join(CURR_DIRECTORY_PATH, os.pardir)), 'datafiles')
TARGET_PATHWAYS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_pathways.csv')
TARGET_ORTHOLOGS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_orthologs.csv')
TARGET_ORGANISMS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_organisms.csv')
TARGET_GENES_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'in', 'target_genes.csv')
TEMP_GENES_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'temp', 'temp_genes.csv')

PATHWAYS_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'pathways_results.csv')
ORTHOLOGS_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'orthologs_results.csv')
GENES_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'genes_results.csv')
TEMP_UNIPROT_IDS_RESULTS_FILEPATH = os.path.join(CURR_DIRECTORY_PATH, 'out', 'temp', 'temp_uniprot_ids_results.csv')

UNIPROT_IDS_FILEPATH = os.path.join(DATAFILES_FILEPATH, 'orthologs_uniprot_ids.csv')
PATHWAYS_UNIPROTS_FILEPATH = os.path.join(DATAFILES_FILEPATH, 'uniprot_ids.csv')
KO_UNIPROTS_FILEPATH = os.path.join(DATAFILES_FILEPATH, 'ko_uniprots.csv')


# Step 1: get a list of the correct pathways for humans
# Step 2: get all KEGG ortholog entries associated with each pathway
# Step 3: get all genes associated with each ortholog
# Step 4: get UniProt IDs for each gene

# returns a dictionary of uniprot_id's grouped by KEGG Ortholog ID's (e.g. {'K03841': 'A0A1U7QCS9', ...})
def find_ortholog_uniprots_by_pathway(target_pathway, target_organisms_filepath):
    print('in kegg_orthologs.find_ortholog_uniprots...')

    # get orthologs for the given pathway
    kegg_orthologs_ids = get_orthologs_by_pathway(target_pathway)
    print(f'Found {len(kegg_orthologs_ids)} KO codes for pathway {target_pathway}')
    print(f'KO codes: {kegg_orthologs_ids}')
    pathway_orthologs_df = pd.DataFrame(columns=['ko_id', 'gene_name']).assign(ko_id=kegg_orthologs_ids)   

    # creates a temporary save file for the genes (so you don't have to start completely over if it crashes)
    if os.path.exists(TEMP_GENES_FILEPATH):
        temp_gene_df = pd.read_csv(TEMP_GENES_FILEPATH)
        unprocessed_orthologs = pathway_orthologs_df.loc[pathway_orthologs_df['ko_id'].isin(temp_gene_df['ko_id']) == False] # only iterate over genes that haven't been already processed

    else:
        temp_gene_df = pd.DataFrame(columns=['ko_id', 'gene_name'])
        unprocessed_orthologs = pathway_orthologs_df

    # get genes for each KEGG ortholog
    # gene_list = []
    for index, row in unprocessed_orthologs.iterrows():
            genes = get_genes_by_ortholog(row["ko_id"]) # gets an array of the gene names
            # gene_list.append(genes)

            # writes the genes to the temp save file
            temp_gene_loop_df = pd.DataFrame([{'ko_id': row["ko_id"], 'gene_name': gene} for gene in genes], columns=['ko_id', 'gene_name'])
            temp_gene_df = pd.concat([temp_gene_df, temp_gene_loop_df])
            temp_gene_df.to_csv(TEMP_GENES_FILEPATH, index=False)


    #write genes to a file
    # ortholog_genes_df = pathway_orthologs_df.assign(gene_name=gene_list).explode('gene_name').dropna(subset=["gene_name"]).drop_duplicates() # cleans up data
    ortholog_genes_df = temp_gene_df.dropna(subset=["gene_name"]).drop_duplicates() # cleans up data
    ortholog_genes_df.to_csv(GENES_RESULTS_FILEPATH, index=False)


    # filter the list of genes for the appropriate organisms (if applicable)
    organisms_df = pd.read_csv(target_organisms_filepath) 
    ortholog_genes_df['kegg_organism_code'] = ortholog_genes_df.apply(get_organism_code_from_row, axis=1)
    filtered_genes = ortholog_genes_df.loc[ortholog_genes_df['kegg_organism_code'].isin(organisms_df['kegg_organism_code'])]

    # creates a temporary save file for the uniprots (so you don't have to start completely over if it crashes)
    if os.path.exists(TEMP_UNIPROT_IDS_RESULTS_FILEPATH):
        temp_uniprot_df = pd.read_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)
        unprocessed_genes = filtered_genes.loc[filtered_genes['gene_name'].isin(temp_uniprot_df['gene_name']) == False] # only iterate over genes that haven't been already processed
    else:
        temp_uniprot_df = pd.DataFrame(columns=['ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
        unprocessed_genes = filtered_genes
    
    # get uniprot ids for each KEGG ortholog
    # id_list = []
    for index, row in unprocessed_genes.iterrows():
        uniprot_id = get_uniprot_ids_by_gene(row["gene_name"])
        # id_list.append(uniprot_id)

        # writes the uniprots to the temp save file
        temp_uniprot_loop_df = pd.DataFrame([{'ko_id': row["ko_id"], 'gene_name': row["gene_name"], 'kegg_organism_code': row["kegg_organism_code"], 'uniprot_id': uniprot_id}], columns=['ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
        temp_uniprot_df = pd.concat([temp_uniprot_df, temp_uniprot_loop_df])
        temp_uniprot_df.to_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH, index=False)
    
    # adds the uniprot_id's to the gene df and filters out genes without uniprot ids
    uniprot_ids_df = temp_uniprot_df.dropna(subset=["uniprot_id"]) 
    # uniprot_ids_df = filtered_genes.assign(uniprot_id=id_list).dropna(subset=["uniprot_id"])       

    #write uniprot ids to a file
    # uniprot_ids_df.to_csv(UNIPROT_IDS_FILEPATH)
    uniprots_filepath = os.path.join(DATAFILES_FILEPATH, "ortholog_uniprots", "pathways", f'uniprot_ids_{target_pathway}.csv')
    uniprot_ids_df.to_csv(uniprots_filepath)

    # condenses the dataframe into a dictionary of ko_id's and their associated uniprot_id's
    grouped_results = (uniprot_ids_df.groupby("ko_id")).agg(uniprot_ids = pd.NamedAgg(column = 'uniprot_id', aggfunc = list)).reset_index()

    results_dictionary = dict(zip(grouped_results['ko_id'], grouped_results['uniprot_ids'])) 

    # removes temp save files
    if os.path.exists(TEMP_UNIPROT_IDS_RESULTS_FILEPATH):
        os.remove(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)
    if os.path.exists(TEMP_GENES_FILEPATH):
        os.remove(TEMP_GENES_FILEPATH)

    return results_dictionary

# returns a list of uniprot_id values
def find_ortholog_uniprots_by_ko_id(target_organisms_filepath, ko_id):
    print('in kegg_orthologs.find_ortholog_uniprots_by_ko_id...')

    # get genes for the specified KEGG ortholog
    genes = get_genes_by_ortholog(ko_id) # gets gene names for the given ko_id (returns list<string>)
    genes_dict = [{'ko_id': ko_id, 'gene_name': gene_name} for gene_name in genes] # TODO: integrate intpo line below after testing
    ortholog_genes_df = pd.DataFrame(genes_dict, columns=['ko_id', 'gene_name'])
    # ortholog_genes_df = pathway_orthologs_df.assign(gene_name=gene_list).explode('gene_name').dropna(subset=["gene_name"]).drop_duplicates()

    # filter the list of genes for the appropriate organisms (if applicable)
    organisms_df = pd.read_csv(target_organisms_filepath) 
    ortholog_genes_df['kegg_organism_code'] = ortholog_genes_df.apply(get_organism_code_from_row, axis=1)
    filtered_genes = ortholog_genes_df.loc[ortholog_genes_df['kegg_organism_code'].isin(organisms_df['kegg_organism_code'])]

    # get uniprot ids for each KEGG ortholog

    # TODO: change this to act like a tempfile
    
    if USING_TEMP_UNIPROT_DATA: 
        temp_uniprot_df = pd.read_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)
        filtered_genes = filtered_genes.loc[filtered_genes['gene_name'].isin(temp_uniprot_df['gene_name']) == False]
    else:
        temp_uniprot_df = pd.DataFrame(columns=[ 'pathway_name', 'pathway_description', 'ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
        
    id_list = []
    for index, row in filtered_genes.iterrows():
        uniprot_id = get_uniprot_ids_by_gene(row["gene_name"])
        id_list.append(uniprot_id)

        if uniprot_id != None:
            temp_uniprot_loop_df = pd.DataFrame([{'pathway_name': row["pathway_name"], 'pathway_description': row["pathway_description"], 'ko_id': row["ko_id"], 'gene_name': row["gene_name"], 'kegg_organism_code': row["kegg_organism_code"], 'uniprot_id': uniprot_id}], columns=[ 'pathway_name', 'pathway_description', 'ko_id', 'gene_name', 'kegg_organism_code', 'uniprot_id'])
            temp_uniprot_df = pd.concat([temp_uniprot_df, temp_uniprot_loop_df])
            temp_uniprot_df.to_csv(TEMP_UNIPROT_IDS_RESULTS_FILEPATH, index=False)
    ####


    uniprot_ids_df = filtered_genes.assign(uniprot_id=id_list).dropna(subset=["uniprot_id"])       

    #write uniprot ids to a file
    # uniprot_ids_df.to_csv(UNIPROT_IDS_FILEPATH)
    uniprots_filepath = os.path.join(DATAFILES_FILEPATH, "ortholog_uniprots", "ko_ids", f'ko_uniprots_{ko_id}.csv')
    uniprot_ids_df.to_csv(uniprots_filepath)

    if os.path.exists(TEMP_UNIPROT_IDS_RESULTS_FILEPATH):
        os.remove(TEMP_UNIPROT_IDS_RESULTS_FILEPATH)

    uniprot_id_array = uniprot_ids_df["uniprot_id"].array
    return uniprot_id_array

# retrieves a list of metabolic pathways matching the PATHWAY_KEYWORDS
# def get_pathways_from_keywords(keyword_array):
#     print('in kegg_orthologs.get_pathways...')

#     # Gets a list of all human pathways
#     human_pathways = REST.kegg_list("pathway", "hsa").read().rstrip().split('\n')

#     target_pathways = []
#     print(f"{len(human_pathways)} pathways found.")

#     # Filter human pathways for target keywords
#     for line in human_pathways:
#         entry, description = line.split("\t")
#         if any(keyword.lower() in description.lower() for keyword in keyword_array):
#             target_pathways.append([entry, description])

#     print(f"{len(target_pathways)} pathways matched the keywords.")

#     # converts the pathway names and descriptions into a Pandas Dataframe
#     target_pathways_df = pd.DataFrame(target_pathways, columns=['pathway_name', 'pathway_description'])

#     return target_pathways_df


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
                # print("Extracted KO code:", ko_code)
            # else:
                # print("KO code not found in the line.")

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

    print(f'gene: {gene_name} --> uniprot_id: {uniprot_id}') # TODO: remove after testing
    return uniprot_id

def get_organism_code_from_row(row):
    return row['gene_name'].split(':')[0]

def read_organisms_from_file(filepath):
    pd.read_csv(filepath) 