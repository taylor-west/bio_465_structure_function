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
    postURL = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"

    email = "alexwalbom@gmail.com"
    title = "practiceSubmission"

    sequences = readDirectoryContents()

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