# If you need to process multiple files, you could use Biopython to parse a PDB structure.
from Bio.PDB import PDBParser
import requests
import json
from urllib.request import urlretrieve
import os
import pandas as pd
import ast

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


def make_expected_cluster_lists_and_find_actual_clusters(invariant_res_df: pd.DataFrame, distance_threshold: float):
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

def filter_interesting_clusters(clusters_df: pd.DataFrame, sequence_separation_threshold: int):
    separated_cluster_col = [None for i in range(len(clusters_df))]
    clusters_df['separated_clusters'] = separated_cluster_col
    for i, row in clusters_df.iterrows():
        pos_list = row['clusters']
        if isinstance(pos_list, str):
            seq_pos = row['seq_pos']
            interesting_clusters_list = []
            pos_list = ast.literal_eval(pos_list)
            for pos in pos_list:
                if abs(seq_pos - pos) > sequence_separation_threshold:
                    interesting_clusters_list.append(pos)

            if len(interesting_clusters_list) > 0:
                clusters_df.at[i, 'separated_clusters'] = str(interesting_clusters_list)
    clusters_df.dropna(inplace=True, ignore_index=True)
    clusters_df.to_csv('../datafiles/cluster_data/interesting_clusters.csv')
    return clusters_df


def find_common_clusters(res_clusters_df: pd.DataFrame):
    invariant_MSA_positions = set(res_clusters_df['MSA_pos'])
    MSA_cluster_positions_list = []
    for i, row in res_clusters_df.iterrows():
        if row['separated_clusters'] is None:
            MSA_cluster_positions_list.append(None)
            continue
        MSA_cluster_positions = []
        separated_cluster = row['separated_clusters']
        if isinstance(separated_cluster, str):
            separated_cluster = ast.literal_eval(separated_cluster)
        for pos in separated_cluster:
            MSA_iteration = 0
            seq_iteration = 0
            seq = row['MSA_sequence']
            for i in range(len(seq)):
                if seq[i] == '-':
                    MSA_iteration += 1
                else:
                    MSA_iteration += 1
                    seq_iteration += 1
                if seq_iteration == pos:
                    MSA_cluster_positions.append(MSA_iteration)
        MSA_cluster_positions_list.append(MSA_cluster_positions)
    res_clusters_df['MSA_cluster_positions'] = MSA_cluster_positions_list
    res_clusters_df.to_csv('../datafiles/cluster_data/filtered_dataframe.csv')
    return res_clusters_df

    # for MSA_pos in invariant_MSA_positions:
    #     same_pos_df = res_clusters_df[res_clusters_df['MSA_pos'] == MSA_pos]





def get_clusters_dataframe():
    invariant_res_df = pd.read_csv('../datafiles/muscle_data/MSA_results.csv')
    invariant_res_df.drop(columns='Unnamed: 0', inplace=True)
    result = make_expected_cluster_lists_and_find_actual_clusters(invariant_res_df, 6)
    # print(result)
    clusters_df = pd.read_csv('../datafiles/cluster_data/clusters.csv')
    result2 = filter_interesting_clusters(clusters_df, 20)
    # print(result2)
    interesting_dataframe = pd.read_csv('../datafiles/cluster_data/interesting_clusters.csv')
    find_common_clusters(interesting_dataframe)

# get_clusters_dataframe()

