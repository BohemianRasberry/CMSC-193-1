def count_fasta_chars(input_file):
  """
  Counts the total number of characters in a FASTA file (excluding whitespaces and newlines).

  Args:
      input_file: Path to the FASTA file.

  Returns:
      The total number of characters in the FASTA file (excluding whitespaces and newlines).
  """
  total_chars = 0
  with open(input_file, 'r') as f:
    for line in f:
      # Skip header lines and empty lines
      if line.startswith('>') or not line.strip():
        continue

      # Count non-whitespace characters
      total_chars += len(line.strip())

  return total_chars

# Example usage
input_file = "tester.sixteenth.fa"
num_chars = count_fasta_chars(input_file)

print(f"Total characters (excluding whitespaces and newlines) in {input_file}: {num_chars}")
