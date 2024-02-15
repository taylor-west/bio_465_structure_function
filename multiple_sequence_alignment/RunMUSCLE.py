
import requests

uniprot_ids = ['P37869', 'A5FUW3', 'A0A411SXA6', 'A0A0P7YND4', 'B6W7V0', 'Q815K8', 'B7GTK2', 'A0A0A7K4W5', 'A0A411DND6', 'C0B602', 'G0EYL2', 'A0A080NRS1', 'A0A1H1MXD9', 'C9RQM6',
               'A0Q5J9', 'A0A0P8B8C6', 'Q5ZTX1', 'A0A8B4R093', 'D5ES25', 'Q88MF9', 'Q0S4I1', 'Q6N5U6', 'A0A0P7VUY9', 'A7AZX6', 'Q08YA3', 'F2R6D4', 'Q31QJ8']

def get_uniprot_entry(entry_id):
    # either of the below URLs will work
    # url = f'https://rest.uniprot.org/uniprotkb/search?query={entry_id}&format=fasta'
    url = f'https://rest.uniprot.org/uniprotkb/{entry_id}.fasta'

    response = requests.get(url)

    if response.status_code == 200:
        return response.text
    else:
        print(f'Error: Unable to retrieve UniProt entry {entry_id}')
        return None

fasta_str = ''
num_successes = 0

for entry_id in uniprot_ids:
    entry_data = get_uniprot_entry(entry_id)

    if entry_data:
        # print(entry_data)
        fasta_str += entry_data
        num_successes += 1

with open('/Users/olivia/Documents/capstone/bio_465_structure_function/multiple_sequence_alignment/output.fasta', 'w') as myFile:
    myFile.write(fasta_str)

    myFile.close()

print('done')
print(f'completed {num_successes} entries')

# Define the endpoint URL
url = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"

# Define the required parameters
email = "osession@byu.edu"
title = "Testing"

# Create a dictionary with the parameters
data = {
    'email': email,
    'title': title,
    'sequence': fasta_str,
}

# Make the POST request
job_response = requests.post(url, data=data)

# Print the response
if job_response.status_code == 200:
    job_id = job_response.text
    print(f'here is the job_id: {job_id}')
    status_response = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{job_id}')
    if status_response.status_code == 200:
        print(status_response.text)
        while status_response.text == 'RUNNING' or status_response.text == 'QUEUED':
            status_response = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{job_id}')
            if status_response.status_code == 200:
                if status_response.text == 'RUNNING' or status_response.text == 'QUEUED':
                    continue
                else:
                    break
            else:
                print('something went wrong')
        alignment_response = requests.get(f'https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{job_id}/out')
        if alignment_response.status_code == 200:
            print(alignment_response.text)
        else:
            print('unable to retrieve final sequence alignment results')
else:
    print('job failed')
    print(job_response.json)


