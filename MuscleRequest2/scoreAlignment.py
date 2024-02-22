import os

from Bio import AlignIO, motifs


#takes only the alignment.alignment object
def make_counter_dictionary_from_alignment_ids(alignment):
    protein_position_counter = {}
    for record in alignment:
        protein_position_counter[record.id] = 0
    return protein_position_counter

def make_position_dictionary_from_alignment_ids(alignment):
    protein_position_counter = {}
    for record in alignment:
        protein_position_counter[record.id] = []
    return protein_position_counter


if __name__ == "__main__":

    # Parse the Clustal format alignment
    alignment = AlignIO.read("alignment.aln", "clustal")

    # Create a Motif object
    # motif = motifs.create(alignment)

    # Get the consensus sequence
    # this is wrong
    # consensus = motif.consensus

    # Define a threshold for conservation score
    threshold = 1.0  # Adjust as needed
    column_annotations = alignment.column_annotations['clustal_consensus']
    # Identify invariant locations
    counter_dictionary = make_counter_dictionary_from_alignment_ids(alignment)
    invariant_locations_dict = make_position_dictionary_from_alignment_ids(alignment)
    alignment_length = alignment.get_alignment_length()

    # Iterate through each column
    # for column_index in range(alignment_length):
    #     # Extract the column (a list of characters for each sequence)
    #     column = alignment[:, column_index]

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

    print(invariant_locations_dict)
    num_invariant = list(column_annotations).count('*')

    for record in alignment:
        alignment_to_uniprot = {}
        uniprot_id = record.id.split("_")[0]
        alignment_to_uniprot[record.id] = uniprot_id




    if os.path.exists(os.path.join(os.getcwd(), "indices")):
        os.unlink(os.path.join(os.getcwd(), "indices"))

    # # Output the list of invariant locations
    # with open("indices", 'a') as outF:
    #     print("Invariant locations (position):")
    #     outF.write("Invariant locations (index):\n")
    #     for position in invariant_locations:
    #         print(position + 1)  # Positions are 0-indexed in Python
    #         outF.write(str(position) + "\n")