# import requests
# import json
import pandas
from Bio.KEGG import REST

PATHWAY_KEYWORDS = ["metabolic", "metabolism"]

# Step 1: get a list of the correct pathways for humans
# Step 2: find all genes/entries associated with those pathways
# Step 3: get KEGG ortholog entry from gene
# Step 4: get all genes associated with that ortholog, making sure to track organisms


# retrieves a list of metabolic pathways matching the PATHWAY_KEYWORDS
def get_pathways():
    print('in kegg.get_pathways...')

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
    target_pathways_df = pandas.DataFrame(target_pathways, columns =['pathway_name', 'pathway_description'])
    pathway_names = target_pathways_df.loc[:, "pathway_name"]
    for pathway_name in pathway_names:
        genes = get_genes_by_pathway(pathway_name)

    return

def get_genes_by_pathway(pathway_name):
    # https://www.biostars.org/p/6224/#6355
    print(f'in kegg.get_genes_by_pathway for pathway_name={pathway_name}...')
    
    # see below link for some information on Bio.KEGG.REST
    # https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html
    # kegg_get() 
    return


def get_genes_by_ortholog():
    print('in kegg.get_genes_by_ortholog...')

    return

def get_ortholog_from_gene():
    print('in kegg.get_ortholog_from_gene...')
    
    return

if __name__ == "__main__":
    get_pathways()



# Exaclt what we want to do, but it is in R instead of Python
# https://support.bioconductor.org/p/109871/