import time
import tracemalloc

def naiveBuildSA(t):
    satups = sorted([(t[i:], i) for i in range(len(t))])
    return list(map(lambda x: x[1], satups))

def lcp(x, y):
    for i in range(min(len(x), len(y))):
        if x[i] != y[i]:
            return i
    return min(len(x), len(y))

def build_generalized_suffix_array(collection, reference):
    combined_string = collection + '$' + reference + '#'
    combined_suffix_array = naiveBuildSA(combined_string)
    
    # Determine which suffixes belong to the collection or the reference
    suffix_type = ['C' if i < len(collection) else 'R' for i in combined_suffix_array]
    
    return combined_suffix_array, suffix_type

def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

# Measure time and space complexity
start_time = time.time()
tracemalloc.start()

# Read sequences from .fa files
collection_sequence = read_fasta('SRR835775_1.clean.fa')
reference_sequence = read_fasta('reference.fa')

# Build the generalized suffix array
gsa, suffix_type = build_generalized_suffix_array(collection_sequence, reference_sequence)

# Stop measuring time and space complexity
end_time = time.time()
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()

# Output the generalized suffix array
'''print("Generalized Suffix Array:")
for index, stype in zip(gsa, suffix_type):
    print(f"{index} ({stype})")'''

# Print time and space complexity
print(f"\nTime taken: {end_time - start_time:.5f} seconds")
print(f"Current memory usage: {current / 1024 / 1024:.4f} MB")
print(f"Peak memory usage: {peak / 1024 / 1024:.4f} MB")
