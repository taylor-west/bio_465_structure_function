import requests
import json
import pandas as pd
import os
from Bio.PDB import PDBParser
from CalculateResidueDistanceWithAverageCoords import download_alphafold_data

uniprot_ids = ['P37869', 'A5FUW3', 'A0A411SXA6', 'A0A0P7YND4', 'B6W7V0', 'Q815K8', 'B7GTK2', 'A0A0A7K4W5', 'A0A411DND6', 'C0B602', 'G0EYL2', 'A0A080NRS1', 'A0A1H1MXD9', 'C9RQM6',
               'A0Q5J9', 'A0A0P8B8C6', 'Q5ZTX1', 'A0A8B4R093', 'D5ES25', 'Q88MF9', 'Q0S4I1', 'Q6N5U6', 'A0A0P7VUY9', 'A7AZX6', 'Q08YA3', 'F2R6D4', 'Q31QJ8']

# def get_uniprot_entry(entry_id):
#     url = f'https://rest.uniprot.org/uniprotkb/{entry_id}.json'
#
#     response = requests.get(url)
#
#     if response.status_code == 200:
#         return json.loads(response.text)
#     else:
#         print(f'Error: Unable to retrieve UniProt entry {entry_id}')
#         return None
#
#
# uniprot_df = pd.DataFrame(columns=['uniprot_id', 'organism_name', 'amino_acid_sequence'])
#
# for entry_id in uniprot_ids:
#     entry_data = get_uniprot_entry(entry_id)
#
#     if entry_data:
#         AAS = entry_data['sequence']['value']
#         organism = entry_data['organism']['scientificName']
#
#         df_row = pd.DataFrame({'uniprot_id': [entry_id], 'organism_name': [organism], 'uniprot_sequence': [AAS]})
#         uniprot_df = pd.concat([uniprot_df, df_row], ignore_index=True)
#
# print(uniprot_df)
#
# residue_codes = pd.read_csv('residue_codes.csv')
#
# pdb_sequences = []
# for i in range(len(uniprot_df)):
#     uniprot_id = uniprot_df.loc[i]['uniprot_id']
#     pdb_filepath = download_alphafold_data(uniprot_id)
#     if pdb_filepath:
#         print(f'download successful for {uniprot_id}')
#         parser = PDBParser()
#
#         # read structure from file
#         structure = parser.get_structure(f'{uniprot_id}', pdb_filepath)
#
#         model = structure[0]
#         chain = model['A']
#         sequence = ''
#         print(chain)
#         residues = chain.get_residues()
#         for residue in residues:
#             row = residue_codes[residue_codes['Three Letter Code'] == residue.get_resname()]
#             single_letter_code = row['Single Letter Code']
#             sequence += list(single_letter_code)[0]
#         print(sequence)
#         pdb_sequences.append(sequence)
#
#     os.remove(pdb_filepath)
#
# uniprot_df['pdb_sequence'] = pdb_sequences
# uniprot_df.to_csv('compare_sequences.csv')

seq_df = pd.read_csv('compare_sequences.csv')

seq_matches = []
for i in range(len(seq_df)):
    uniprot_seq = seq_df.loc[i]['uniprot_sequence']
    pdb_seq = seq_df.loc[i]['pdb_sequence']
    if uniprot_seq == pdb_seq:
        seq_matches.append(True)
        print(True)
    else:
        seq_matches.append(False)
        print(False)

seq_df['sequences_match'] = seq_matches


