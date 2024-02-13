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
    # either of the below URLs will work
    # url = f'https://rest.uniprot.org/uniprotkb/search?query={entry_id}&format=fasta'
    # url = f'https://rest.uniprot.org/uniprotkb/{entry_id}.fasta'
    url = f'https://www.ebi.ac.uk/proteins/api/proteins/{entry_id}'

    response = requests.get(url)

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

if __name__ == "__main__":
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

    with open(os.path.join(os.getcwd, "uniprotEntries"), 'w') as myFile:
        myFile.write(sequences)

        myFile.close(sequences)

    data = {
        "email": email,
        "title": title,
        "sequence": sequences
    }

    responseCode = requests.post(postURL, data=data)
    jobID = responseCode.text

    getURL1 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/resulttypes/{jobID}"

    status = checkStatus(jobID)

    # resultType = requests.get(getURL1)
    # print(resultType)

    format = "aln-clustalw"

    getURL2 = f"https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{jobID}/{format}"
    alignment = requests.get(getURL2)
    print(alignment.text)
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