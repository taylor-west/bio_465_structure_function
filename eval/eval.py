import pandas as pd
import os
import json
import requests
from urllib.request import urlretrieve

from Bio.PDB import PDBParser

from cluster_analysis_3d.CalculateResidueDistanceWithDataframeInput import find_clusters, filter_interesting_clusters, find_common_clusters

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

        # is_key = []
        # for pos in group['seq_pos']:
        #     if pos in all_key_positions:
        #         is_key.append(True)
        #     else:
        #         is_key.append(False)
        # is_key = pd.Series(is_key)
        key_positions_list = [all_key_positions] * len(only_key_rows_df)
        # Add a new column 'key_positions' to the filtered DataFrame
        only_key_rows_df['key_positions'] = key_positions_list
        aggregated_df = pd.concat([aggregated_df, only_key_rows_df], ignore_index=True)
    return aggregated_df

def make_expected_cluster_lists_for_all_invariant_residues(invariant_res_df: pd.DataFrame, distance_threshold: float):

        residue_codes = pd.read_csv('../datafiles/cluster_data/residue_codes.csv')

        filepath = '../datafiles/pdb_files/'
        os.makedirs(filepath, exist_ok=True)

        uniprot_ids = set(invariant_res_df['uniprot_id'])
        cluster_col = [None for i in range(len(invariant_res_df))]
        invariant_res_df['clusters'] = cluster_col
        invariant_res_df['expected_residues'] = None
        # going through every uniprot id in our dataframe
        grouped_by_uniprot = invariant_res_df.groupby('uniprot_id')
        for uniprot_id, group in grouped_by_uniprot:
            filtered_by_id_df = invariant_res_df[invariant_res_df['uniprot_id'] == uniprot_id]
            entry_id = filtered_by_id_df['entry_id'].iloc[0]
            uniprot_filename = entry_id + ".txt"
            key_positions = read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, uniprot_filename)
            # Store key_positions in a list
            act_positions = key_positions.get('ACT_SITE', [])
            bind_positions = key_positions.get('BINDING', [])
            all_key_positions = act_positions + bind_positions

            clusters_dict = {}
            url = f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}'

            response = requests.get(url)

            alphafold_data = None

            if response.status_code == 200:
                alphafold_data = response.text
            else:
                print(f'Error: Unable to retrieve AlphaFold entry for {uniprot_id}')

            if alphafold_data:
                # get the pdb url to download the file
                alphafold_json = json.loads(alphafold_data)
                pdb_url = alphafold_json[0]['pdbUrl']
                print(pdb_url)

                filename = f'alphafold_{uniprot_id}.pdb'
                urlretrieve(pdb_url, filepath + filename)

                if os.path.exists(filepath + filename):

                    # create parser
                    parser = PDBParser()

                    # read structure from file
                    # TODO: update the id in this next line
                    structure = parser.get_structure('ENOLASE (A0A2P1UMH5)', filepath + filename)

                    # set up empty dictionary to track clusters
                    locs_3d_dict = {}
                    model = structure[0]
                    chain = model['A']

                    # create a dictionary to hold the 3d positions of every uniprot residue
                    uniprot_residue_3d_positions = {}
                    for position in all_key_positions:
                        residue = chain[position]
                        uniprot_residue_3d_positions[position] = residue['CA']


                    for i, row in filtered_by_id_df.iterrows():
                        residue = chain[row['seq_pos']]
                        csv_map_row = residue_codes[residue_codes['Three Letter Code'] == residue.get_resname()]
                        single_letter_code = list(csv_map_row['Single Letter Code'])[0]

                        if single_letter_code == row['residue']:
                            locs_3d_dict[row['seq_pos']] = residue['CA']

                    for pos in locs_3d_dict:
                        clusters_list = []
                        expected_residues_list = []
                        # calculate distance between current invariant residue and all key residues for this protein
                        for uniprot_residue, position_3d in uniprot_residue_3d_positions.items():
                            distance = locs_3d_dict[pos] - position_3d
                            if distance <= distance_threshold:
                                if distance != 0:
                                    expected_residues_list.append(uniprot_residue)

                        for pos2 in locs_3d_dict:
                            if pos != pos2:
                                # calculate distance between residues
                                distance = locs_3d_dict[pos] - locs_3d_dict[pos2]
                                if distance <= distance_threshold:
                                    clusters_list.append(pos2)

                        if len(clusters_list) > 0:
                            index = invariant_res_df.loc[(invariant_res_df['uniprot_id'] == uniprot_id) & (
                                        invariant_res_df['seq_pos'] == pos), 'clusters'].index[0]
                            invariant_res_df.at[index, 'clusters'] = clusters_list
                            invariant_res_df.at[index, 'expected_residues'] = expected_residues_list
            for file in os.listdir(filepath):
                os.remove(f'{filepath}{file}')
        invariant_res_df = invariant_res_df.drop(columns=['MSA_sequence'])
        invariant_res_df.to_csv('../datafiles/cluster_data/clusters.csv', index=False)
        return invariant_res_df


def analyze_key_locations(dataframe):
    key_results_df = get_key_rows_only(dataframe)
    key_results_df['common_positions'] = None
    key_results_df['num_key_present'] = 0
    key_results_df['percent_key_present'] = None
    for index, row in key_results_df.iterrows():
        cluster_positions = set(row['clusters'])
        key_positions = set(row['key_positions'])
        common_positions = list(cluster_positions.intersection(key_positions))
        key_results_df.at[index, 'common_positions'] = common_positions
        num_key_present = len(common_positions)
        percent_key_present = len(common_positions) / len(cluster_positions)
        key_results_df.at[index, 'num_key_present'] = num_key_present
        key_results_df.at[index, 'percent_key_present'] = percent_key_present
    return key_results_df

def evaluate_results(dataframe):
    key_df = get_key_rows_only(dataframe)
    print(key_df)
    results_df = analyze_key_locations(key_df)
    # Modifies the original DataFrame
    results_df.drop('MSA_sequence', axis=1, inplace=True)
    results_df.to_csv(os.path.join(PATH_TO_EVAL_FILES, 'analyzed_clusters.csv'))

invariant_res_df = pd.read_csv('../datafiles/muscle_data/MSA_results.csv')
# clusters_df = pd.read_csv('../datafiles/cluster_data/clusters.csv')
# invariant_res_df.drop(columns='Unnamed: 0', inplace=True)

# result = find_clusters(invariant_res_df, 10)

#result2 = cluster.filter_interesting_clusters(clusters_df, 20)
# print(result.dtypes)
# evaluate_results(result)

make_expected_cluster_lists_for_all_invariant_residues(invariant_res_df, distance_threshold=7)







