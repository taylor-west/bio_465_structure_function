import requests
import os
    
def get_uniprot_file(entry_id):
    # gets the fasta file for a protein in UniProt by id
    url = f'https://rest.uniprot.org/uniprotkb/{entry_id}.fasta'

    response = requests.get(url)

    if response.status_code == 200:
        return response.text
    else:
        return None
    
file_name = "TEXT-FILE-WITH-PROTEIN-IDS" # replace with file containing protein ids from kegg pathways
prot_ids = []
for line in open(file_name,'r').readlines():
    prot_ids += [line.strip()]

for id in prot_ids:
    # make folder for ortholog files
    directory = f'Ortholog_Sequences/{id}'
    if not os.path.exists(directory):
        os.mkdir(directory)

    # opens file resulting from oma_api_calls.py
    f = open(f'Orthologs/{id}_orthologs.tsv','r')

    # orthologs read from the last column of aforementioned file
    orthologs = []
    index = 0
    for line in f.readlines():
        if index != 0:
            orthologs += [line.strip().split('\t')[2]]
        index += 1

    f.close()

    for ortholog in orthologs:
        entry_data = get_uniprot_file(ortholog)
        if entry_data:
            f = open(f'{directory}/{ortholog}.fasta','w')
            f.write(entry_data)
            f.close()
        else:
            # creates file for proteins that did not have a valid UniProt id listed in OMA
            f = open(f'unavailable_proteins.tsv','a')
            f.write(f'{ortholog}\n')
        