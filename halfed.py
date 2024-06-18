def cut_fasta_half_lines(input_file, output_file):
  """
  Cuts a FASTA file in half (number of lines) and writes it to a new file.

  Args:
      input_file: Path to the input FASTA file.
      output_file: Path to the output FASTA file containing half the lines.
  """
  with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    # Count total lines (excluding empty lines)
    total_lines = sum(1 for line in f_in if line.strip())

    # Handle empty file
    if total_lines == 0:
      return

    # Determine half the number of lines (rounded down)
    half_lines = total_lines // 2

    # Reset file pointer to the beginning
    f_in.seek(0)

    # Write the first half of the lines (excluding empty lines)
    written_lines = 0
    while written_lines < half_lines:
      line = f_in.readline()
      if line.strip():
        f_out.write(line)
        written_lines += 1

# Example usage
input_file = "tester.eighth.fa"
output_file = "tester.sixteenth.fa"
cut_fasta_half_lines(input_file, output_file)

print(f"First half (up to lines) from {input_file} written to {output_file}")