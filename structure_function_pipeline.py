import os
import ortholog_identification

HelpOutput = """
# python kegg_orthologs.py -p file.txt -m file.txt -n file.txt

Required Parameters:
-k  <file>  This is a csv file containing a list of keywords that will be used to filter pathways in KEGG.

-p  <file>  This is a string naming the KEGG pathway ID that will be used to retrieve orthologs in KEGG.

-o   <file>  This is a csv file containing the KEGG orthologs corresponding to the genes whose uniprot_id's will be retrieved.
-g   <file>  This is a csv file of genes for which to get the uniprot_id's.

Example: 
"""

if __name__ == "__main__":
    # TODO: for testing only. Read variables in from command line arguements later.
    kegg_pathway = '2hsa00051'
    target_organisms_filepath = os.path.join(os.path.cwd(), 'target_organisms.csv')
    target_ko_id = 'K03841'
    ##


    # generate the list of uniprot ids for the given pathway and organisms (and optionally a ko_id value)
    # returns a list of uniprot_id's in an array
    if target_ko_id == None:
        # takes a string naming the KEGG Pathway ID value and a filepath pointing to a csv file containing the target organisms
        uniprots = ortholog_identification.find_ortholog_uniprots(kegg_pathway, target_organisms_filepath)
    else:
        # takes a filepath pointing to a csv file containing the target organisms and the KEGG Ortholog ID value for the genes of interest
        uniprots = ortholog_identification.find_ortholog_uniprots_by_ko_id(target_organisms_filepath, target_ko_id)

    print(f'located {len(uniprots)} UniProt IDs')