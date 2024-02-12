
# If you need to process multiple files, you could use Biopython to parse a PDB structure.

from Bio.PDB import PDBParser
import requests
import json
from urllib.request import urlretrieve
import os

uniprot_id = 'A0A2P1UMH5'

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

    # TODO: update this filepath for the proper place on your machine, we will have to set this up to be a parameter in the future
    filepath = '/Users/olivia/PycharmProjects/PDBparsing/'
    filename = f'alphafold_{uniprot_id}.pdb'
    urlretrieve(pdb_url, filepath + filename)

    if os.path.exists(filepath + filename):

        # create parser
        parser = PDBParser()

        # read structure from file
        structure = parser.get_structure('ENOLASE (A0A2P1UMH5)', filepath + filename)
        #'AF-A0A2P1UMH5-F1-model_v4.pdb'

        desired_residue = 'MET'
        conserved_positions = {1: None, 37: None, 87: None}

        model = structure[0]
        chain = model['A']

        residue_pos = 1
        print(len(chain))
        for residue in chain:
            if residue.get_resname() == desired_residue and residue_pos in conserved_positions.keys():
                conserved_positions[residue_pos] = residue['CA']
            residue_pos += 1

        print(conserved_positions)

        close_residues = {}
        # this dict will map which residues are actually close to each other
        # if all the lists remain empty by the end, of the double for loop,
        # then we can conclude that none of the conserved residues are clustered together
        # else we can conclude that some of the residues are clustered together and we can
        # see exactly which ones were in close proximity in 3d space
        for position in conserved_positions:
            close_residues[position] = []

        max_distance = 2 # later we will turn this into a parameter that can be adjusted to see how it affects our results

        for pos in conserved_positions:
            for pos2 in conserved_positions:
                if pos != pos2:
                    # calculate distance between residues
                    distance = conserved_positions[pos] - conserved_positions[pos2]
                    if distance <= max_distance:
                        close_residues[pos].append(pos2)

        print(close_residues)


# this code comes from a biostars post

# # this example uses only the first residue of a single chain.
# # it is easy to extend this to multiple chains and residues.
# for residue1 in chain:
#     for residue2 in chain:
#         if residue1 != residue2:
#             # compute distance between CA atoms
#             try:
#                 distance = residue1['CA'] - residue2['CA']
#             except KeyError:
#                 ## no CA atom, e.g. for H_NAG
#                 continue
#             if distance < 6:
#                 print(residue1, residue2, distance)
#         # stop after first residue
#         break