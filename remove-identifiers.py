def remove_odd_lines(input_filename, output_filename):
  """
  Removes every odd line from an input FASTA file and writes the result to an output file.

  Args:
      input_filename (str): Path to the input FASTA file.
      output_filename (str): Path to the output FASTA file where even lines will be written.
  """

  with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
    for line_number, line in enumerate(input_file, start=1):
      if line_number % 2 == 0:  # Check for even line number
        output_file.write(line)

# Example usage
input_filename = "tester.fa"
output_filename = "tester.clean.fa"
remove_odd_lines(input_filename, output_filename)
print(f"Even lines written to: {output_filename}")
