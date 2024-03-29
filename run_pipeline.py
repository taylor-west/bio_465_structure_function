import os
import sys
from ortholog_identification.kegg_orthologs import find_ortholog_uniprots_by_ko_id
from ortholog_identification.filter_orthologs_with_binding_sites import find_uniprot_gene_collections
from multiple_sequence_alignment.muscle import multiple_sequence_alignment
from cluster_analysis_3d import CalculateResidueDistanceWithDataframeInput
from eval import eval
import pandas as pd

CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "datafiles")

if __name__ == "__main__":
    num_cmdline_args = len(sys.argv)-1
    target_organisms_filepath = sys.argv[1]
    target_ko_id = sys.argv[2]

    if not os.path.exists(PATH_TO_DATAFILES):
        os.makedirs(PATH_TO_DATAFILES)

    # takes a filepath pointing to a csv file containing the target organisms and the KEGG Ortholog ID value for the genes of interest
    # returns a list of uniprot ids for the given KEGG Ortholog ID value
    uniprots = find_ortholog_uniprots_by_ko_id(target_organisms_filepath, target_ko_id)
    kegg_results_df = find_uniprot_gene_collections(f'datafiles/ortholog_uniprots/ko_uniprots_{target_ko_id}.csv', f'datafiles/ortholog_uniprots/temp.tsv')
    print(kegg_results_df)
    unique_orthologs = set(kegg_results_df['ko_id'])
    for ortholog in unique_orthologs:
        subset_df = kegg_results_df[kegg_results_df['ko_id'] == ortholog]
        uniprots = list(subset_df['uniprot_id'])
        print(f'located {len(uniprots)} UniProt IDs for ortholog {ortholog}')
        print(uniprots)
        # generates a multiple sequence alignment from the uniprot ids
        msa_df = multiple_sequence_alignment(uniprots, ortholog)
        print("done with alignment")
        invariant_res_df = CalculateResidueDistanceWithDataframeInput.make_expected_cluster_lists_and_find_actual_clusters(ortholog, 7)
        print(f"done with distance calculations for {ortholog}")
        # eval.evaluate_results(invariant_res_df) # TODO: decide if needed
        print("done with evaluation")
