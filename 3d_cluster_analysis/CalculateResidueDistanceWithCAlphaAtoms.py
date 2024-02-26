
# If you need to process multiple files, you could use Biopython to parse a PDB structure.

from Bio.PDB import PDBParser
import requests
import json
from urllib.request import urlretrieve
import os
import pandas as pd


with open('../datafiles/3d-cluster-test-data.txt', 'r') as file:
    dict_str = file.read()
    invariant_locs_dict = eval(dict_str)


def find_clusters(invariant_locs_dict: dict, distance_threshold: float):
    residue_codes = pd.read_csv('residue_codes.csv')

    clusters_for_prots = {}
    clusters_dict = {}
    filepath = '../datafiles/pdb_files/'
    for uniprot_id, locs in invariant_locs_dict.items():
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
                # TODO: update this next line's id
                structure = parser.get_structure('ENOLASE (A0A2P1UMH5)', filepath + filename)
                #'AF-A0A2P1UMH5-F1-model_v4.pdb'

                # set up empty dictionary to track clusters
                locs_3d_dict = {}
                model = structure[0]
                chain = model['A']

                for invariant_res in locs:
                    clusters_dict[invariant_res[0]] = []
                    print(invariant_res)
                    residue = chain[invariant_res[0]]
                    row = residue_codes[residue_codes['Three Letter Code'] == residue.get_resname()]
                    single_letter_code = list(row['Single Letter Code'])[0]
                    print(single_letter_code)
                    if single_letter_code == invariant_res[1]:
                        locs_3d_dict[invariant_res[0]] = residue['CA']


                print(locs_3d_dict)
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
                                clusters_dict[pos].append(pos2)

        clusters_for_prots[uniprot_id] = clusters_dict
        for file in os.listdir(filepath):
            os.remove(f'{filepath}{file}')

    return clusters_for_prots

result = find_clusters(invariant_locs_dict, 5)
print(result)