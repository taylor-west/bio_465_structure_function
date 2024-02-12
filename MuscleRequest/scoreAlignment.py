from Bio import AlignIO
from Bio.Align import AlignInfo

# Parse the Clustal format alignment
alignment = AlignIO.read("alignment.clustal", "clustal")

summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()

# Define a threshold for conservation score
threshold = 1.0  # Adjust as needed

# Identify invariant locations
invariant_locations = []
for i, char in enumerate(consensus):
    if char != '-':
        column = alignment[:, i]
        if column.count(char) / len(column) >= threshold:
            invariant_locations.append(i)

# Output the list of invariant locations
print("Invariant locations (position):")
for position in invariant_locations:
    print(position + 1)  # Positions are 0-indexed in Python