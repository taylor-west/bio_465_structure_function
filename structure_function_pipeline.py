import sys

import multiple_sequence_alignment.muscle
import ortholog_identification.kegg_orthologs

if __name__ == "__main__":
    # uses KOID and organisms or all three not sure
    # TODO update this statement with somethinf for the command line args (if statement):
    KOID = sys.argv[1]
    organisms = sys.argv[2]
    
    # generate the list of uniprot ids for the given pathway
    uniprot_ids = ortholog_identification.kegg_orthologs.get_uniprot_ids_for_pathways(KOID, organisms)
    msa_df = multiple_sequence_alignment.muscle.multiple_sequence_alignment()

    # the other option is only use a pathway
    pathway = sys.argv[1]
    kegg_orthologs = ortholog_identification.kegg_orthologs.get_uniprot_ids_for_pathways(pathway)
    for ortholog in kegg_orthologs:
        # extract the uniprot ids for all of that ortholog
        extracted_ids = []
        msa_df = multiple_sequence_alignment.muscle.multiple_sequence_alignment(extracted_ids)
        #continue to 3d
