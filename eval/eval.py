import pandas as pd
import os
import json
import requests
from urllib.request import urlretrieve

from Bio.PDB import PDBParser

import cluster_analysis_3d.CalculateResidueDistanceWithDataframeInput

CWD = os.getcwd()
PATH_TO_DATAFILES = os.path.join(CWD, "../datafiles")
PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")


def read_uniprot_file_to_analyze_active_sites(directory, filename):
    important_positions = {}
    act_str = 'ACT_SITE'
    bind_str = 'BINDING'
    important_positions[act_str] = []
    important_positions[bind_str] = []
    for file in os.listdir(directory):
        if file == filename:
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as inF:
                for line in inF:
                    line = line.strip()
                    if line.startswith('FT'):
                        entry = line.split()
                        if entry[1] == act_str:
                            important_positions[act_str].append(int(entry[2]))
                        if entry[1] == bind_str:
                            important_positions[bind_str].append(int(entry[2]))
    return important_positions

# returns a dataframe that only includes entries with uniprot residues as the seq_pos
def get_key_rows_only(dataframe):
    grouped_df = dataframe.groupby('entry_id')
    aggregated_df = pd.DataFrame()
    for entry_id, group in grouped_df:
        entry_id = group['entry_id'].iloc[0]
        filename = entry_id + ".txt"
        key_positions = read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, filename)
        # Store key_positions in a list
        act_positions = key_positions.get('ACT_SITE', [])
        bind_positions = key_positions.get('BINDING', [])
        all_key_positions = act_positions + bind_positions
        # makes a dataframe that only includes the rows with key positions
        only_key_rows_df = group[group.seq_pos.isin(all_key_positions)].copy()

        key_positions_list = [all_key_positions] * len(only_key_rows_df)
        # Add a new column 'key_positions' to the filtered DataFrame
        only_key_rows_df['key_positions'] = key_positions_list
        aggregated_df = pd.concat([aggregated_df, only_key_rows_df], ignore_index=True)
    return aggregated_df


def analyze_key_locations(dataframe):
    # key_results_df = get_key_rows_only(dataframe)
    key_results_df = dataframe
    key_results_df['common_positions'] = None
    key_results_df['num_key_present'] = 0
    key_results_df['percent_of_cluster_expected'] = None
    key_results_df['percent_expected_found'] = None
    for index, row in key_results_df.iterrows():
        cluster_positions = set(row['clusters'])
        expected_positions = set(row['expected_residues'])
        common_positions = list(cluster_positions.intersection(expected_positions))
        if not common_positions:
            key_results_df.at[index, 'common_positions'] = []
        key_results_df.at[index, 'common_positions'] = common_positions
        num_key_present = len(common_positions)
        percent_of_cluster_expected = len(common_positions) / len(cluster_positions)
        if len(expected_positions) == 0:
            percent_of_cluster_expected = -1
            percent_expected_found = -1
        else:
            percent_expected_found = len(common_positions) / len(expected_positions)
        key_results_df.at[index, 'num_key_present'] = num_key_present
        key_results_df.at[index, 'percent_of_cluster_expected'] = percent_of_cluster_expected
        key_results_df.at[index, 'percent_expected_found'] = percent_expected_found
    return key_results_df

def evaluate_results(dataframe):
    # key_df = get_key_rows_only(dataframe)
    # print(key_df)
    results_df = analyze_key_locations(dataframe)
    # Modifies the original DataFrame
    # results_df.drop('MSA_sequence', axis=1, inplace=True)
    results_df.to_csv(os.path.join(PATH_TO_EVAL_FILES, 'analyzed_clusters.csv'))

invariant_res_df = pd.read_csv('../datafiles/muscle_data/MSA_results.csv')
# clusters_df = pd.read_csv('../datafiles/cluster_data/clusters.csv')
# invariant_res_df.drop(columns='Unnamed: 0', inplace=True)

result = cluster_analysis_3d.CalculateResidueDistanceWithDataframeInput.make_expected_cluster_lists_and_find_actual_clusters(invariant_res_df, 10)

#result2 = cluster.filter_interesting_clusters(clusters_df, 20)
# print(result.dtypes)
evaluate_results(result)







