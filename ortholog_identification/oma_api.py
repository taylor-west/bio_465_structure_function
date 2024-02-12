# requires pip install omadb
from omadb import Client
import os 

c = Client()

# will eventually take protein ids from kegg pathways
file_name = "TEXT-FILE-WITH-PROTEIN-IDS" # replace with file containing protein ids from kegg pathways
prot_ids = []
for line in open(file_name,'r').readlines():
    prot_ids += [line.strip()]

# make folder for ortholog files
directory = 'Orthologs/'
if not os.path.exists(directory):
    os.mkdir(directory)

for id in prot_ids:
    r = c.proteins[id]  # Can also be called as c.proteins.info(prot_id)

    # opens output file
    f = open(f'Orthologs/{id}_orthologs.tsv','w')

    # inserts column names
    f.write('Species\tOriginal Gene\tOrtholog\n')

    # includes species code, id of starting protein, and ids of all the starting proteins' orthologs
    for ortholog in r.orthologs:
        f.write(f'{ortholog["species"]["code"]}\t{id}\t{ortholog["canonicalid"]}\n')