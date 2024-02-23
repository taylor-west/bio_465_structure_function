import os
from Bio import AlignIO
#ALL OF THIS RUNS IN main.py BUT WE ARE NOT SURE HOW WE ARE SEPARATING THE FILES YET

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

# returns a dictionary of uniprotIDs and their correspinding sequences
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

if __name__ == "__main__":
    sequenceDictionary = read_uniprot_files(os.path.join(os.getcwd(), "uniprotEntries"))
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