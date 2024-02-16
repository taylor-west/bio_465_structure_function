import os

from Bio import AlignIO, motifs

# Parse the Clustal format alignment
alignment = AlignIO.read("alignment.aln", "clustal")

# Create a Motif object
motif = motifs.create(alignment)

# Get the consensus sequence
consensus = motif.consensus

# Define a threshold for conservation score
threshold = 1.0  # Adjust as needed

# Identify invariant locations
invariant_locations = []
for i, char in enumerate(consensus):
    if char != '-':
        column = alignment[:, i]
        if column.count(char) / len(column) >= threshold:
            invariant_locations.append(i)

if os.path.exists(os.path.join(os.getcwd(), "indices")):
    os.unlink(os.path.join(os.getcwd(), "indices"))

# Output the list of invariant locations
with open("indices", 'a') as outF:
    print("Invariant locations (position):")
    outF.write("Invariant locations (index):\n")
    for position in invariant_locations:
        print(position + 1)  # Positions are 0-indexed in Python
        outF.write(str(position) + "\n")