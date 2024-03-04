import os
import sys
from ortholog_identification.kegg_orthologs import find_ortholog_uniprots_by_pathway, find_ortholog_uniprots_by_ko_id
import multiple_sequence_alignment.muscle

if __name__ == "__main__":
    num_cmdline_args = len(sys.argv)-1
    kegg_pathway = sys.argv[1]
    # kegg_pathway = '2hsa00051' # TODO: remove once commandline args work
    target_organisms_filepath = sys.argv[2]
    # target_organisms_filepath = os.path.join(os.path.cwd(), 'target_organisms.csv') # TODO: remove once commandline args work
    target_ko_id = (sys.argv[3] if num_cmdline_args > 3 else None)
    # target_ko_id = 'K03841' # TODO: remove once commandline args work


    # generate the list of uniprot ids for the given pathway and organisms (and optionally a ko_id value)
    if target_ko_id == None:
        # takes a string naming the KEGG Pathway ID value and a filepath pointing to a csv file containing the target organisms
        # returns a dictionary of uniprot_id's grouped by KEGG Ortholog ID's (e.g. {'K03841': 'A0A1U7QCS9', ...})
        ortholog_uniprot_dict = find_ortholog_uniprots_by_pathway(kegg_pathway, target_organisms_filepath)
        print(f'located {len(ortholog_uniprot_dict)} orthologs with {len(ortholog_uniprot_dict.values())} UniProt IDs')

        # loops through each KEGG Ortholog and generates a multiple sequence alignment for each
        for ortholog in ortholog_uniprot_dict:
            # extract the uniprot ids for all of that ortholog
            extracted_ids = ortholog_uniprot_dict[ortholog]
            msa_df = multiple_sequence_alignment.muscle.multiple_sequence_alignment(extracted_ids)
            #continue to 3d
    else:
        # takes a filepath pointing to a csv file containing the target organisms and the KEGG Ortholog ID value for the genes of interest
        # returns a list of uniprot ids for the given KEGG Ortholog ID value
        uniprots = find_ortholog_uniprots_by_ko_id(target_organisms_filepath, target_ko_id)
        print(f'located {len(uniprots)} UniProt IDs')

        # generates a multiple sequence alignment from the uniprot ids
        msa_df = multiple_sequence_alignment.muscle.multiple_sequence_alignment()

    #continue to 3d