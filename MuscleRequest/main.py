import requests
import time
import os
from Bio import AlignIO

# makes a dictionary to keep track of what position we are at for each protein
def make_counter_dictionary_from_alignment_ids(alignment):
    protein_position_counter = {}
    for record in alignment:
        protein_position_counter[record.id] = 0
    return protein_position_counter

# makes a dict to keep track of the locations of every invariant residue in each protein
def make_position_dictionary_from_alignment_ids(alignment):
    protein_position_counter = {}
    for record in alignment:
        protein_position_counter[record.id] = []
    return protein_position_counter

# retrieves the status of Muscle job
def getStatus(jobID):
    statusURL = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{jobID}"
    return requests.get(statusURL)

# checks if Muscle job is finished periodically
def checkStatus(jobID):
    keepChecking = True
    while keepChecking:
        time.sleep(5)
        status = getStatus(jobID)
        if status.text == "FINISHED":
            keepChecking = False
        else:
            print("Waiting for results...")

# given a set of accession values, retrieves uniprot flatfiles for each protein
def get_uniprot_entry(entry_id):
    url = f'https://www.ebi.ac.uk/proteins/api/proteins/{entry_id}'
    response = requests.get(url, headers={"Accept": "text/x-flatfile"})
    if response.status_code == 200:
        return response.text
    else:
        print(f'Error: Unable to retrieve UniProt entry {entry_id}')
        return None

def readDirectoryContents(folderPath):
    fileList = os.listdir(folder_path)
    sequences = ""
    for fileName in fileList:
        filePath = os.path.join(folder_path, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            sequences += inF.read()
            sequences += "\n"
    return sequences

def remove_files_in_subfolder(folder_path):
    # Get the list of files in the subfolder
    files_to_remove = os.listdir(folder_path)
    # if os.path.exists(os.path.join(os.getcwd(), "alignment.aln")):
    #     os.unlink(os.path.join(os.getcwd(), "alignment.aln"))
    # Iterate through the files and remove them
    for file_name in files_to_remove:
        file_path = os.path.join(folder_path, file_name)
        try:
            os.unlink(file_path)
        except Exception as e:
            print(f"Error removing {file_path}: {e}")


# returns a dictionary of uniprotIDs and their corresponding sequences
def read_uniprot_files(directory):
    uniprot_data = {}
    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as file:
                uniprot_id = None
                sequence = ''
                for line in file:
                    line = line.strip()
                    if line.startswith('ID'):
                        uniprot_id = line.split()[1]
                    elif line.startswith('SQ'):
                        sequence_lines = []
                        for line in file:
                            if line.startswith('//'):
                                break
                            sequence_lines.append(line.strip())
                        sequence = ''.join(sequence_lines).replace(' ', '')
                        break
                if uniprot_id and sequence:
                    uniprot_data[uniprot_id] = sequence
    return uniprot_data

def retrieve_uniprot_id(filepath):
    with open(filepath, 'r') as file:
        uniprot_id = None
        for line in file:
            line = line.strip()
            if line.startswith('AC'):
                uniprot_id = line.split()[1][0:-1]
                break
    return uniprot_id


if __name__ == "__main__":
    folder_path = os.path.join(os.getcwd(), "../uniprotEntries")
    if len(os.listdir(folder_path)) > 0:
        remove_files_in_subfolder(folder_path)
    uniprot_ids = ['P34949', 'A5A6K3', 'G3RFM0', 'G7MYC5', 'A0A2K6DHS4', 'A0A096NMS2', 'A0A2K5JTJ0', 'U3CWX3', 'A0A2K6T9B3',
        'A0A1U7UAV0', 'A0A8B7E9X4', 'Q924M7', 'A0A6P5Q4Y5', 'Q68FX1', 'G3I837', 'A0A1U7QCS9', 'A0A8C6R367', 'G5AL67',
        'A0A8B7VRU0', 'A0A1S3ETZ6']

    for entry_id in uniprot_ids:
        entry_data = get_uniprot_entry(entry_id)
        if entry_data:
            first_line = entry_data.split()
            file_name = first_line[1]
            filePath = os.path.join(os.getcwd(), "../uniprotEntries", f"{file_name}.txt")
            with open(filePath, 'w') as inF:
                inF.write(entry_data)

    postURL = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"

    email = "alexwalbom@gmail.com"
    title = "practiceSubmission"

    sequences = readDirectoryContents(folder_path)

    with open(os.path.join(os.getcwd(), "sequenceFile.txt"), 'w') as myFile:
        myFile.write(sequences)

    data = {
        "email": email,
        "title": title,
        "sequence": sequences
    }

    responseCode = requests.post(postURL, data=data)
    jobID = responseCode.text

    getURL1 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/resulttypes/{jobID}"

    # returns when the job is complete
    checkStatus(jobID)

    format = "aln-clustalw"

    getURL2 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{jobID}/{format}"
    alignment = requests.get(getURL2)
    file = open(os.path.join(os.getcwd(), "alignment.aln"), "w")
    file.write(alignment.text)
    file.close()

    sequenceDictionary = read_uniprot_files(os.path.join(os.getcwd(), "../uniprotEntries"))
    # Parse the Clustal format alignment
    alignment = AlignIO.read("alignment.aln", "clustal")

    # Define a threshold for conservation score
    threshold = 1.0  # Adjust as needed
    column_annotations = alignment.column_annotations['clustal_consensus']
    # Identify invariant locations
    counter_dictionary = make_counter_dictionary_from_alignment_ids(alignment)
    invariant_locations_dict = make_position_dictionary_from_alignment_ids(alignment)
    alignment_length = alignment.get_alignment_length()

    for i, char in enumerate(column_annotations):
        meets_threshold = False
        column = alignment[:, i]
        most_frequent_element = max(set(column), key=column.count)
        frequency = column.count(most_frequent_element) / len(column)
        if frequency >= threshold:
            meets_threshold = True
        # update the index
        for seq in alignment:
            if seq[i] != '-':
                counter_dictionary[seq.id] += 1
            if meets_threshold:
                invariant_locations_dict[seq.id].append((counter_dictionary[seq.id], seq[i]))

    uniprot_invariant_locs_dict = {}
    for key, val in invariant_locations_dict.items():
        uniprot_id = retrieve_uniprot_id(f'{os.path.join(os.getcwd(), "../uniprotEntries")}/{key}.txt')
        uniprot_invariant_locs_dict[uniprot_id] = val

    print(uniprot_invariant_locs_dict)
    # this dict is what will be passed on to the 3d_cluster_analysis