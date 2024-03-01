
# If you need to process multiple files, you could use Biopython to parse a PDB structure.

from Bio.PDB import PDBParser
import requests
import json
from urllib.request import urlretrieve
import os
import pandas as pd
import ast



def find_clusters(invariant_res_df: pd.DataFrame, distance_threshold: float):
    residue_codes = pd.read_csv('residue_codes.csv')

    clusters_for_prots = {}
    filepath = '../datafiles/pdb_files/'

    uniprot_ids = set(invariant_res_df['uniprot_id'])
    cluster_col = [None for i in range(len(invariant_res_df))]
    invariant_res_df['clusters'] = cluster_col
    for uniprot_id in uniprot_ids:
        filtered_by_id_df = invariant_res_df[invariant_res_df['uniprot_id'] == uniprot_id]
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

                for i, row in filtered_by_id_df.iterrows():
                    residue = chain[row['seq_pos']]
                    csv_map_row = residue_codes[residue_codes['Three Letter Code'] == residue.get_resname()]
                    single_letter_code = list(csv_map_row['Single Letter Code'])[0]
                    if single_letter_code == row['residue']:
                        locs_3d_dict[row['seq_pos']] = residue['CA']

                # this dict will map which residues are actually close to each other
                # if all the lists remain empty by the end, of the double for loop,
                # then we can conclude that none of the conserved residues are clustered together
                # else we can conclude that some of the residues are clustered together and we can
                # see exactly which ones were in close proximity in 3d space

                for pos in locs_3d_dict:
                    clusters_list = []
                    for pos2 in locs_3d_dict:
                        if pos != pos2:
                            # calculate distance between residues
                            distance = locs_3d_dict[pos] - locs_3d_dict[pos2]
                            if distance <= distance_threshold:
                                clusters_list.append(pos2)
                    if len(clusters_list) > 0:
                        index = invariant_res_df.loc[(invariant_res_df['uniprot_id'] == uniprot_id) & (invariant_res_df['seq_pos'] == pos), 'clusters'].index[0]
                        invariant_res_df.at[index, 'clusters'] = clusters_list

        clusters_for_prots[uniprot_id] = clusters_dict
        for file in os.listdir(filepath):
            os.remove(f'{filepath}{file}')
    invariant_res_df = invariant_res_df.drop(columns='Unnamed: 0')
    invariant_res_df.to_csv('../datafiles/clusters.csv')
    return clusters_for_prots

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
    # clean up the dataframe
    clusters_df.drop(columns=['Unnamed: 0'], inplace=True)
    clusters_df.dropna(inplace=True, ignore_index=True)
    clusters_df.to_csv('interesting_clusters.csv')
    return clusters_df

def find_common_clusters(res_clusters_df: pd.DataFrame):
    already_searched = []
    for i, row in res_clusters_df.iterrows():
        MSA_pos = row['MSA_pos']
        if MSA_pos in already_searched:
            continue
        already_searched.append(MSA_pos)
        same_pos_df = res_clusters_df[res_clusters_df['MSA_pos'] == MSA_pos]
    pass

invariant_res_df = pd.read_csv('../datafiles/MSA_results.csv')
clusters_df = pd.read_csv('../datafiles/clusters.csv')
result = find_clusters(invariant_res_df, 10)
result2 = filter_interesting_clusters(clusters_df, 20)
print(result)
print(result2)