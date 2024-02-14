import requests
import json
import pandas as pd
import os
from urllib.request import urlretrieve
from Bio.PDB import PDBParser

uniprot_ids = ['P37869', 'A5FUW3', 'A0A411SXA6', 'A0A0P7YND4', 'B6W7V0', 'Q815K8', 'B7GTK2', 'A0A0A7K4W5', 'A0A411DND6',
               'C0B602', 'G0EYL2', 'A0A080NRS1', 'A0A1H1MXD9', 'C9RQM6',
               'A0Q5J9', 'A0A0P8B8C6', 'Q5ZTX1', 'A0A8B4R093', 'D5ES25', 'Q88MF9', 'Q0S4I1', 'Q6N5U6', 'A0A0P7VUY9',
               'A7AZX6', 'Q08YA3', 'F2R6D4', 'Q31QJ8']

def get_uniprot_entry(entry_id):
    url = f'https://rest.uniprot.org/uniprotkb/{entry_id}.json'

    response = requests.get(url)

    if response.status_code == 200:
        return json.loads(response.text)
    else:
        print(f'Error: Unable to retrieve UniProt entry {entry_id}')
        return None

def get_uniprot_data_df(entry_ids: list):
    uniprot_df = pd.DataFrame(columns=['uniprot_id', 'organism_name', 'amino_acid_sequence'])

    for entry_id in uniprot_ids:
        entry_data = get_uniprot_entry(entry_id)

        if entry_data:
            AAS = entry_data['sequence']['value']
            organism = entry_data['organism']['scientificName']

            df_row = pd.DataFrame({'uniprot_id': [entry_id], 'organism_name': [organism], 'uniprot_sequence': [AAS]})
            uniprot_df = pd.concat([uniprot_df, df_row], ignore_index=True)

    return uniprot_df

def get_pdb_ids(uniprot_id, organism_name):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.structure_determination_methodology",
                        "operator": "exact_match",
                        "value": "experimental"
                    }
                },
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "full_text",
                            "parameters": {
                                "value": uniprot_id
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "full_text",
                            "parameters": {
                                "value": organism_name
                            }
                        }
                    ],
                    "label": "full_text"
                }
            ],
            "logical_operator": "and"
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 25
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }
    response = requests.post(url, json=query)
    if response.status_code == 200:
        data = response.json()
        print(data)
        pdb_ids = [hit['identifier'] for hit in data['result_set']]
        return pdb_ids
    else:
        print(response.status_code)
        print(f"Error: Unable to retrieve PDB IDs for UniProt ID {uniprot_id}")
        return []


def download_pdb_files(pdb_ids: list, uniprot_id: str):
    successful_filepaths = []
    base_url = "https://files.rcsb.org/download/"
    for pdb_id in pdb_ids:
        url = f"{base_url}{pdb_id}.pdb"
        filepath = f'../3d_cluster_analysis/{uniprot_id}_pdb_files/'
        filename = f'protein_data_bank_{pdb_id}.pdb'
        if os.path.exists(filepath):
            urlretrieve(url, filepath + filename)
            successful_filepaths.append(f'{filepath}{filename}')
        else:
            os.mkdir(f'../3d_cluster_analysis/{uniprot_id}_pdb_files/')
            urlretrieve(url, filepath + filename)
            successful_filepaths.append(f'{filepath}{filename}')
    return successful_filepaths

def find_pdb_sequences(pdb_filepaths, uniprot_id):
    residue_codes = pd.read_csv('residue_codes.csv')

    pdb_sequences = []
    for pdb_filepath in pdb_filepaths:
        parser = PDBParser()

        # read structure from file
        structure = parser.get_structure(f'{uniprot_id}', pdb_filepath)

        model = structure[0]
        chain = model['A']
        sequence = ''
        residues = chain.get_residues()
        for residue in residues:
            row = residue_codes[residue_codes['Three Letter Code'] == residue.get_resname()]
            single_letter_code = row['Single Letter Code']
            sequence += list(single_letter_code)[0]
        print(sequence)
        pdb_sequences.append(sequence)

        os.remove(pdb_filepath)
    return pdb_sequences

# Example usage
uniprot_df = get_uniprot_data_df(uniprot_ids)
pdb_ids = get_pdb_ids(uniprot_df.loc[0]['uniprot_id'], uniprot_df.loc[0]['organism_name'])
if pdb_ids:
    filepaths = download_pdb_files(pdb_ids, uniprot_df.loc[0]['uniprot_id'])
    pdb_sequences = find_pdb_sequences(filepaths, uniprot_df.loc[0]['uniprot_id'])

    for pdb_seq in pdb_sequences:
        if uniprot_df.loc[0]['uniprot_sequence'] == pdb_seq:
            print(True)



