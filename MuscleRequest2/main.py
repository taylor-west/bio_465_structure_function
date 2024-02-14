import requests
import time
import os

def getStatus(jobID):
    statusURL = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{jobID}"
    return requests.get(statusURL)

def checkStatus(jobID):
    keepChecking = True
    while keepChecking:
        time.sleep(5)
        status = getStatus(jobID)
        if status.text == "FINISHED":
            keepChecking = False
            
def get_uniprot_entry(entry_id):
    url = f'https://www.ebi.ac.uk/proteins/api/proteins/{entry_id}'

    response = requests.get(url, headers={ "Accept" : "text/x-flatfile"})

    if response.status_code == 200:
        return response.text
    else:
        print(f'Error: Unable to retrieve UniProt entry {entry_id}')
        return None

def readDirectoryContents():
    currentDirectory = os.getcwd()
    folderName = "uniprotEntries"
    subdirectoryPath = os.path.join(currentDirectory, folderName)
    fileList = os.listdir(subdirectoryPath)
    sequences = ""

    for fileName in fileList:
        filePath = os.path.join(subdirectoryPath, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            sequences += inF.read()
            sequences += "\n"

    return sequences

def remove_files_in_subfolder(folder_path):
    # Get the list of files in the subfolder
    files_to_remove = os.listdir(folder_path)
    if os.path.exists(os.path.join(os.getcwd(), "alignment.clustal")):
        os.unlink(os.path.join(os.getcwd(), "alignment.clustal"))

    # Iterate through the files and remove them
    for file_name in files_to_remove:
        file_path = os.path.join(folder_path, file_name)
        try:
            os.unlink(file_path)
        except Exception as e:
            print(f"Error removing {file_path}: {e}")

if __name__ == "__main__":
    folder_path = os.path.join(os.getcwd(), "uniprotEntries")
    remove_files_in_subfolder(folder_path)
    uniprot_ids = ['P37869', 'A5FUW3', 'A0A411SXA6', 'A0A0P7YND4', 'B6W7V0', 'Q815K8', 'B7GTK2', 'A0A0A7K4W5', 'A0A411DND6', 'C0B602', 'G0EYL2', 'A0A080NRS1', 'A0A1H1MXD9', 'C9RQM6',
               'A0Q5J9', 'A0A0P8B8C6', 'Q5ZTX1', 'A0A8B4R093', 'D5ES25', 'Q88MF9', 'Q0S4I1', 'Q6N5U6', 'A0A0P7VUY9', 'A7AZX6', 'Q08YA3', 'F2R6D4', 'Q31QJ8']

    fasta_str = ''
    num_successes = 0

    for entry_id in uniprot_ids:
        entry_data = get_uniprot_entry(entry_id)
        if entry_data:
            # print(entry_data)
            currentDirectory = os.getcwd()
            filePath = os.path.join(os.getcwd(), "uniprotEntries", f"{entry_id}.txt")
            with open(filePath, 'w') as inF:
                inF.write(entry_data)
            fasta_str += entry_data
            num_successes += 1

    
    # with open(os.path.join(os.getcwd, "uniprotEntries"), 'w') as myFile:
    #     myFile.write(fasta_str)

    #     myFile.close()

    print('done')
    print(f'completed {num_successes} entries')
    postURL = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"

    email = "alexwalbom@gmail.com"
    title = "practiceSubmission"

    sequences = readDirectoryContents()

    with open(os.path.join(os.getcwd(), "uniprotEntries", "sequenceFile.txt"), 'w') as myFile:
        myFile.write(sequences)

        myFile.close()

    data = {
        "email": email,
        "title": title,
        "sequence": sequences,
        "tree": "tree1"
    }

    responseCode = requests.post(postURL, data=data)
    jobID = responseCode.text

    getURL1 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/resulttypes/{jobID}"

    status = checkStatus(jobID)

    format = "aln-clustalw"

    getURL2 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{jobID}/{format}"
    alignment = requests.get(getURL2)
    f = open("alignment.clustal", "a")
    f.write(alignment.text)
    f.close()




    # Check the status code and handle the response
    if responseCode.status_code == 200:
        print("Request was successful!")
        print(responseCode.text)
    else:
        print(f"Request failed with status code {responseCode.status_code}")
        print(responseCode.text)