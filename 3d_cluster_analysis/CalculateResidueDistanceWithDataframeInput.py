
# If you need to process multiple files, you could use Biopython to parse a PDB structure.

from Bio.PDB import PDBParser
import requests
import json
from urllib.request import urlretrieve
import os
import pandas as pd


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
                    for pos2 in locs_3d_dict:
                        if pos != pos2:
                            # calculate distance between residues
                            distance = locs_3d_dict[pos] - locs_3d_dict[pos2]
                            if distance <= distance_threshold:
                                clusters_list = \
                                filtered_by_id_df.loc[filtered_by_id_df['seq_pos'] == pos, 'clusters'].iloc[0]
                                if clusters_list is not None:
                                    clusters_list.append(pos2)
                                else:
                                    clusters_list = [pos2]
                                # Update the DataFrame directly
                                invariant_res_df.loc[(invariant_res_df['uniprot_id'] == uniprot_id) & (
                                            invariant_res_df['seq_pos'] == pos), 'clusters'] = clusters_list

                # for pos in locs_3d_dict:
                #     for pos2 in locs_3d_dict:
                #         if pos != pos2:
                #             # calculate distance between residues
                #             distance = locs_3d_dict[pos] - locs_3d_dict[pos2]
                #             if distance <= distance_threshold:
                #                 clusters_list = list(filtered_by_id_df[filtered_by_id_df['seq_pos'] == pos]['clusters'])[0]
                #                 if clusters_list is not None:
                #                     clusters_list.append(pos2)
                #                     # filtered_by_id_df[filtered_by_id_df['seq_pos'] == pos]['clusters'] = clusters_list
                #                     invariant_res_df.loc[(invariant_res_df['uniprot_id'] == uniprot_id) & (invariant_res_df['seq_pos'] == pos), 'clusters'] = clusters_list
                #                 else:
                #                     clusters_list = [pos2]
                #                     invariant_res_df.loc[(invariant_res_df['uniprot_id'] == uniprot_id) & (invariant_res_df['seq_pos'] == pos), 'clusters'] = clusters_list
                #                     # filtered_by_id_df[filtered_by_id_df['seq_pos'] == pos]['clusters'] = clusters_list

        clusters_for_prots[uniprot_id] = clusters_dict
        for file in os.listdir(filepath):
            os.remove(f'{filepath}{file}')

    invariant_res_df.to_csv('../datafiles/clusters.csv')
    return clusters_for_prots

def filter_interesting_clusters(uniprot_clusters_dict: dict, sequence_separation_threshold: int):
    interesting_clusters_prots = {}
    for uniprot_id, clusters in uniprot_clusters_dict.items():
        interesting_clusters_dict = {}
        for key_position, value_list in clusters.items():
            for position in value_list:
                if abs(key_position - position) > sequence_separation_threshold:
                    if key_position in interesting_clusters_dict:
                        interesting_clusters_dict[key_position].append(position)
                    else:
                        interesting_clusters_dict[key_position] = []
                        interesting_clusters_dict[key_position].append(position)
        interesting_clusters_prots[uniprot_id] = interesting_clusters_dict
    return interesting_clusters_prots

invariant_res_df = pd.read_csv('../datafiles/MSA_results.csv')
result = find_clusters(invariant_res_df, 10)
result2 = filter_interesting_clusters(result, 20)
print(result)
print(result2)