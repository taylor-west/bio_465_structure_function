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
    # TODO: for testing only. Read variable in from command line arguements later.
    kegg_pathway = '2hsa00051'
    target_organisms_filepath = os.path.join(os.path.cwd(), 'target_organisms.csv')
    ##


    # generate the list of uniprot ids for the given pathway and organisms
    # takes a string naming the KEGG pathway id, and a filepath pointing to a csv file containing the target organisms
    ortholog_identification.get_uniprot_ids_for_pathways(kegg_pathway, target_organisms_filepath)