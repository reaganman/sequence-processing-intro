fasta_file="$1"
samples_dir="$2"
# Get number of samples
num_samples=$(grep -c ">" "$fasta_file")

# Get variable sites
python3 get_variable_sites.py "$fasta_file" "$samples_dir"

