import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def read_fasta(file_path):
    return SeqIO.parse(file_path, "fasta")

def translate_sequence(nucleotide_sequence):
    return str(Seq(nucleotide_sequence).translate())

def make_amino_acid_alignment(nucleotide_sequences):
    amino_acid_sequences = [SeqRecord(Seq(translate_sequence(str(seq.seq))), id=seq.id) for seq in nucleotide_sequences]
    amino_acid_alignment = MultipleSeqAlignment(amino_acid_sequences)
    return amino_acid_alignment

def find_variable_positions(alignment):
    alignment_length = len(alignment[0])
    variable_positions = []

    for position in range(alignment_length):
        column = alignment[:, position]
        if len(set(column)) > 1:
            variable_positions.append(position + 1)  # Adjust for 1-based indexing

    return variable_positions

def create_sample_directories(nucleotide_sequences, output_directory):
    os.makedirs(output_directory, exist_ok=True)

    for seq in nucleotide_sequences:
        sample_directory = os.path.join(output_directory, seq.id)
        os.makedirs(sample_directory, exist_ok=True)

        # Save DNA sequence to file
        dna_sequence_file = os.path.join(sample_directory, "dna_sequence.fasta")
        SeqIO.write(seq, dna_sequence_file, "fasta")

        # Save amino acid sequence to file
        amino_acid_sequence = translate_sequence(str(seq.seq))
        amino_acid_record = SeqRecord(Seq(amino_acid_sequence), id=seq.id)
        amino_acid_sequence_file = os.path.join(sample_directory, "amino_acid_sequence.fasta")
        SeqIO.write(amino_acid_record, amino_acid_sequence_file, "fasta")

def save_translated_sequences(nucleotide_sequences, output_directory):
    # Save all translated sequences to one file
    all_sequences = [SeqRecord(Seq(translate_sequence(str(seq.seq))), id=seq.id) for seq in nucleotide_sequences]
    all_sequences_file = os.path.join(output_directory, "all_translated_sequences.fasta")
    SeqIO.write(all_sequences, all_sequences_file, "fasta")

def save_alignment(alignment):
    alignment_file = "amino_acid_alignment.fasta"
    SeqIO.write(alignment, alignment_file, "fasta")


def get_unique_dates(nucleotide_sequences):
    unique_dates = set(seq.description.split('_')[1] for seq in nucleotide_sequences)
    return list(unique_dates)

def main():
    parser = argparse.ArgumentParser(description="Translate the original alignment and identify variable positions.")
    parser.add_argument("input_file", help="Path to the input alignment file in FASTA format.")
    parser.add_argument("output_directory", help="Path to the output directory.")
    args = parser.parse_args()

    nucleotide_sequences = list(read_fasta(args.input_file))

    # Make nucleotide variable positions
    nucleotide_alignment = MultipleSeqAlignment(nucleotide_sequences)
    variable_positions_nucleotide = find_variable_positions(nucleotide_alignment)

    # Make amino acid alignment and find variable positions
    amino_acid_alignment = make_amino_acid_alignment(nucleotide_sequences)
    variable_positions_amino_acid = find_variable_positions(amino_acid_alignment)

    print("Variable Nucleotide Positions:", variable_positions_nucleotide)
    print("Variable Amino Acid Positions:", variable_positions_amino_acid)

    # Create directories for each sample containing DNA and amino acid sequence files
    create_sample_directories(nucleotide_sequences, args.output_directory)

    # Save all translated sequences to one file
    save_translated_sequences(nucleotide_sequences, args.output_directory)

    # Save amino acid alignment to file
    save_alignment(amino_acid_alignment)

    # Get unique dates
    unique_dates = get_unique_dates(nucleotide_sequences)

    # Save information to log.txt
    with open(os.path.join(args.output_directory, "log.txt"), "w") as log_file:
        log_file.write(f"Number of samples: {len(nucleotide_sequences)}\n")
        log_file.write(f"Number of unique sampling dates: {len(unique_dates)}\n")
        log_file.write(f"Number of variable sites in DNA sequence: {len(variable_positions_nucleotide)}\n")
        log_file.write(f"Number of variable sites in protein sequence: {len(variable_positions_amino_acid)}\n")
        log_file.write("Unique Sampling Dates: {}\n".format(", ".join(unique_dates)))

if __name__ == "__main__":
    main()
