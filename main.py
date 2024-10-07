from Bio import SeqIO, Seq
from Bio.Seq import Seq
import numpy as np

protein_names = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# Generate all possible dicodons
all_dicodons = [protein1 + protein2 for protein1 in protein_names for protein2 in protein_names]

codon_frequencies = []
dicodon_frequencies = []
seq_id = []

fasta_files = [
    "./viruses/data/bacterial1.fasta",
    "./viruses/data/bacterial2.fasta",
    "./viruses/data/bacterial3.fasta",
    "./viruses/data/bacterial4.fasta",
    "./viruses/data/mamalian1.fasta",
    "./viruses/data/mamalian2.fasta",
    "./viruses/data/mamalian3.fasta",
    "./viruses/data/mamalian4.fasta"
]

def find_protein_sequences(sequence):
    protein_sequence = ""
    index = 0
    while True:
        start_index = 0
        stop_index = 0
        index = sequence.find("ATG", index)
        if index != -1:
            start_index = index
            index += 3

            find_result = sequence.find("TAG", index)
            if find_result > 0:
                stop_index = find_result
            find_result = sequence.find("TAA", index)
            if find_result > 0 and find_result < stop_index:
                stop_index = find_result
            find_result = sequence.find("TGA", index)
            if find_result > 0 and find_result < stop_index:
                stop_index = find_result

            if stop_index == 0:
                break
            if (stop_index - start_index + 3) % 3 != 0:
                continue

            if stop_index - start_index < 100:
                seq = Seq(sequence[start_index:stop_index + 3])
                protein_sequence += seq.translate(table="Bacterial", to_stop=True)
        else:
            break
    return protein_sequence


for fasta_file in fasta_files:
    with open(fasta_file, "r") as file:
        codon_frequencies_seq = np.zeros(len(protein_names))
        dicodon_frequencies_seq = np.zeros(len(all_dicodons))
        total_proteins_in_sequence = 0

        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            reverse_complement = record.reverse_complement()
            reverse_complement_sequence = str(reverse_complement.seq)
            seq_id.append(str(record.id))

            # Find protein sequences and update frequencies
            protein_sequence = find_protein_sequences(sequence) + find_protein_sequences(reverse_complement_sequence)

            for protein in protein_sequence:
                if protein in protein_names:
                    protein_index = protein_names.index(protein)
                    codon_frequencies_seq[protein_index] += 1
                    total_proteins_in_sequence += 1

            # Find dicodon sequences and update frequencies
            for start_index in range(len(sequence)):
                if start_index + 6 <= len(sequence):
                    codon1 = sequence[start_index:start_index + 3]
                    codon2 = sequence[start_index + 3:start_index + 6]
                    dicodon = Seq(codon1).translate(table="Bacterial", to_stop=True) + Seq(codon2).translate(
                        table="Bacterial", to_stop=True)
                    if dicodon in all_dicodons:
                        dicodon_index = all_dicodons.index(dicodon)
                        dicodon_frequencies_seq[dicodon_index] += 1

            for start_index in range(len(reverse_complement_sequence)):
                if start_index + 6 <= len(reverse_complement_sequence):
                    codon1 = reverse_complement_sequence[start_index:start_index + 3]
                    codon2 = reverse_complement_sequence[start_index + 3:start_index + 6]
                    dicodon = Seq(codon1).translate(table="Bacterial", to_stop=True) + Seq(codon2).translate(
                        table="Bacterial", to_stop=True)
                    if dicodon in all_dicodons:
                        dicodon_index = all_dicodons.index(dicodon)
                        dicodon_frequencies_seq[dicodon_index] += 1

        # Divide the frequencies by the total protein count in the sequence
        if total_proteins_in_sequence > 0:
            codon_frequencies_seq /= total_proteins_in_sequence
            dicodon_frequencies_seq /= total_proteins_in_sequence

        codon_frequencies.append(codon_frequencies_seq)
        dicodon_frequencies.append(dicodon_frequencies_seq)

# Calculate the distances matrix using Euclidean distance for codons
num_sequences = len(codon_frequencies)
codon_distances_matrix = np.zeros((num_sequences, num_sequences))

for i in range(num_sequences):
    for j in range(i, num_sequences):
        distance = np.linalg.norm(codon_frequencies[i] - codon_frequencies[j])
        codon_distances_matrix[i, j] = distance
        codon_distances_matrix[j, i] = distance

# Calculate the distances matrix using Euclidean distance for dicodons
num_sequences = len(dicodon_frequencies)
dicodon_distances_matrix = np.zeros((num_sequences, num_sequences))

for i in range(num_sequences):
    for j in range(i, num_sequences):
        distance = np.linalg.norm(dicodon_frequencies[i] - dicodon_frequencies[j])
        dicodon_distances_matrix[i, j] = distance
        dicodon_distances_matrix[j, i] = distance

# Format the distances matrices in PHYLIP format
matrix_size_codon = len(codon_distances_matrix)
matrix_size_dicodon = len(dicodon_distances_matrix)
phylip_matrix_codon = []
phylip_matrix_dicodon = []

for i in range(matrix_size_codon):
    row_codon = [f"{codon_distances_matrix[i, j]:.4f}" for j in range(matrix_size_codon)]
    row_codon.insert(0, seq_id[i])
    phylip_matrix_codon.append(" ".join(row_codon))

for i in range(matrix_size_dicodon):
    row_dicodon = [f"{dicodon_distances_matrix[i, j]:.4f}" for j in range(matrix_size_dicodon)]
    row_dicodon.insert(0, seq_id[i])
    phylip_matrix_dicodon.append(" ".join(row_dicodon))

# Output distances matrices
with open("matrix_codon.phylip", "w") as outfile:
    outfile.write(f"{matrix_size_codon}\n")
    outfile.write("\n".join(phylip_matrix_codon))

with open("matrix_dicodon.phylip", "w") as outfile:
    outfile.write(f"{matrix_size_dicodon}\n")
    outfile.write("\n".join(phylip_matrix_dicodon))
