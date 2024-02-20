import requests
import re
import os
import pandas as pd
from Bio.KEGG import REST

READ_PATHWAY_FROM_FILE = False
READ_ORTHOLOGS_FROM_FILE = True
READ_GENES_FROM_FILE = True
USING_TEMP_GENE_DATA = True
PATHWAY_KEYWORDS = ["metabolic", "metabolism"]

CURR_DIRECTORY_PATH = os.path.dirname(__file__)
TARGET_PATHWAYS_FILENAME = 'target_pathways.csv'
PATHWAY_ORTHOLOGS_FILENAME = 'pathway_orthologs.csv'
ORTHOLOG_GENES_FILENAME = 'ortholog_genes.csv'
TEMP_GENE_FILENAME = 'temp_gene_data.csv'
UNIPROT_IDS_FILENAME = 'uniprot_ids.csv'

# Step 1: get a list of the correct pathways for humans
# Step 2: get all KEGG ortholog entries associated with each pathway
# Step 3: get all genes associated with each ortholog
# Step 4: get UniProt IDs for each gene

def get_uniprot_ids_for_pathways():
    print('in kegg_orthologs.get_uniprot_ids_for_pathways...')

    # get pathways data
    if READ_PATHWAY_FROM_FILE:
       pathways_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'working', TARGET_PATHWAYS_FILENAME)) 
       
    else:
        pathways = get_pathways()
        pathways_df = pd.DataFrame(pathways, columns=['pathway_name', 'pathway_description'])

        # write pathways to a file (FOR TESTING/DEV)
        pathways_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', TARGET_PATHWAYS_FILENAME))

    # get orthologs for each pathway
    if READ_ORTHOLOGS_FROM_FILE:
       pathway_orthologs_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'working', PATHWAY_ORTHOLOGS_FILENAME)) 
       
    else:
        kegg_ids_list = []
        for index, row in pathways_df.iterrows():
            kegg_ids = get_orthologs_by_pathway(row['pathway_name'])
            kegg_ids_list.append(kegg_ids)
        pathway_orthologs_df = pathways_df.assign(ko_id=kegg_ids_list).explode('ko_id').dropna(subset=["ko_id"]).drop_duplicates()       

        # write orthologs to a file (FOR TESTING/DEV)
        pathway_orthologs_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'working', PATHWAY_ORTHOLOGS_FILENAME))




    # get genes for each KEGG ortholog
    if READ_GENES_FROM_FILE:
        ortholog_genes_df = pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'working', ORTHOLOG_GENES_FILENAME)) 
    else:
        # ortholog_genes_df = df = pathway_orthologs_df.reindex(columns=["pathway_name", 'pathway_description', "ko_id", "gene_name", "gene_organism"])
        gene_list = []
        temp_gene_df = (pd.DataFrame(columns=['ko_id', 'gene_name']) if USING_TEMP_GENE_DATA else pd.read_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', TEMP_GENE_FILENAME)))
        temp_gene_df
        for index, row in pathway_orthologs_df.iterrows():
            if not (USING_TEMP_GENE_DATA and (row["ko_id"] in temp_gene_df['ko_id'])):
                genes = get_genes_by_ortholog(row["ko_id"], pathway_orthologs_df) # gets an array of the gene names
                gene_list.append(genes)

                # TODO: remove this after testing
                temp_df = pd.DataFrame([{'ko_id': row["ko_id"], 'gene_name': gene} for gene in genes], columns=['ko_id', 'gene_name'])
                temp_gene_df = pd.concat([temp_gene_df, temp_df])
                temp_gene_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'temp', 'working', TEMP_GENE_FILENAME))
                ####

        ortholog_genes_df = pathway_orthologs_df.assign(gene_name=gene_list).explode('gene_name').dropna(subset=["gene_name"])       
        
        #write genes to a file
        ortholog_genes_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', ORTHOLOG_GENES_FILENAME))


    # get uniprot ids for each KEGG ortholog
    # uniprot_ids_df = ortholog_genes_df.reindex(columns = ortholog_genes_df.columns.tolist() + ['uniprot_id'])
    id_list = []
    for index, row in ortholog_genes_df.iterrows():
        uniprot_id = get_uniprot_ids_by_gene(row["gene_name"])
        id_list.append(uniprot_id)
        # row['uniprot_id']=uniprot_id

        # TODO: remove this after testing
        # target_organisms = ['hsa', 'ptr', 'pps', 'ggo', 'pon', 'nle', 'hmh', 'mcc']
        # temp_df = pd.DataFrame([{'ko_id': row["ko_id"], 'gene_name': gene} for gene in genes], columns=['ko_id', 'gene_name'])
        # temp_gene_df = pd.concat([temp_gene_df, temp_df])
        # temp_gene_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'temp', TEMP_GENE_FILENAME))
    uniprot_ids_df = ortholog_genes_df.assign(uniprot_ids=id_list).dropna(subset=["uniprot_ids"])       

    #write uniprot ids to a file
    uniprot_ids_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'working', UNIPROT_IDS_FILENAME))
    
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
def get_genes_by_ortholog(kegg_ortholog_id, temp_save_df):
    print(f'in kegg_orthologs.get_genes_by_ortholog for kegg_ortholog_id={kegg_ortholog_id}...')

    ortholog_genes_response = REST.kegg_find("genes", kegg_ortholog_id).read().rstrip().split('\n')
    
    # # This is another option if we need to include gene description data in the final dataframe. It is significantly slower though.
    # ortholog_genes_response_df = pd.DataFrame([kegg_ortholog_id, [line.split("\t")[0] for line in ortholog_genes_response]], columns=['ko_id', 'gene_name'])
    # temp_save_df.append(ortholog_genes_response_df)
    # temp_save_df.to_csv(os.path.join(CURR_DIRECTORY_PATH, 'data', 'temp', 'temp_gene_data.csv'))
    
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
    # Regular expression pattern to match UniProt ID
    pattern = r'UniProt:\s+(\S+)'
    # Extracting UniProt ID using regular expression
    match = re.search(pattern, gene_response)
    if match:
        uniprot_id = match.group(1)
        print("UniProt ID:", uniprot_id)
    else:
        uniprot_id = None
        print("UniProt ID not found.")
    return uniprot_id

if __name__ == "__main__":
    get_uniprot_ids_for_pathways()