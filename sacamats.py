import time
import tracemalloc

class SuffixArrayData:
    def __init__(self, collection, reference):
        self.collection = collection
        self.reference = reference
        self.combined_string = collection + '$' + reference + '#'
        self.sa = []  # Suffix Array
        self.isa = []  # Inverted Suffix Array
        self.lcp = []  # Longest Common Prefix

def read_fasta(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        return sequence
    except IOError as e:
        print(f"Error reading {filename}: {e}")
        return None

def build_suffix_array(s):
    n = len(s)
    suffix_array = sorted(range(n), key=lambda i: s[i:])
    return suffix_array

def compute_lcp_array(text, sa):
    n = len(text)
    rank = [0] * n
    lcp = [0] * n
    for i, suffix in enumerate(sa):
        rank[suffix] = i
    h = 0
    for i in range(n):
        if rank[i] > 0:
            j = sa[rank[i] - 1]
            while i + h < n and j + h < n and text[i + h] == text[j + h]:
                h += 1
            lcp[rank[i]] = h
            if h > 0:
                h -= 1
    return lcp

def preprocess_reference(reference):
    return reference + '$'

def compute_ecms(collection):
    sa = build_suffix_array(collection)
    isa = [0] * len(collection)
    for i, suffix in enumerate(sa):
        isa[suffix] = i
    lcp = compute_lcp_array(collection, sa)
    return sa, isa, lcp

def create_buckets(lcp, sa):
    buckets = {}
    current_bucket = []
    prev_lcp = 0
    for i, suffix in enumerate(sa):
        if lcp[i] != prev_lcp:
            if current_bucket:
                buckets[prev_lcp] = current_bucket
            current_bucket = [suffix]
            prev_lcp = lcp[i]
        else:
            current_bucket.append(suffix)
    if current_bucket:
        buckets[prev_lcp] = current_bucket
    return buckets

def sort_suffixes_within_buckets(buckets, sa):
    for key in buckets.keys():
        buckets[key].sort(key=lambda x: len(sa[x:]))
    return buckets

def induce_suffixes(buckets, sa):
    sa = induce_suffixes_first_scan(sa, buckets)
    sa = induce_suffixes_second_scan(sa, buckets)
    return sa

def induce_suffixes_first_scan(sa, buckets):
    for key in sorted(buckets.keys()):
        for suffix in buckets[key]:
            if suffix > 0 and sa[suffix - 1] > sa[suffix]:
                sa[sa[suffix - 1]] = suffix - 1
    return sa

def induce_suffixes_second_scan(sa, buckets):
    for key in sorted(buckets.keys(), reverse=True):
        for suffix in buckets[key]:
            if suffix > 0 and sa[suffix - 1] < sa[suffix]:
                sa[sa[suffix - 1]] = suffix - 1
    return sa

def build_generalized_suffix_array(collection_sequence, reference_sequence):
    combined_string = collection_sequence + '#' + reference_sequence + '$'
    sa = build_suffix_array(combined_string)
    isa = [0] * len(combined_string)
    for i, suffix in enumerate(sa):
        isa[suffix] = i
    lcp = compute_lcp_array(combined_string, sa)
    suffix_type = ['C' if idx < len(collection_sequence) else 'R' for idx in sa]
    return sa, isa, lcp, suffix_type

# Measure time and space complexity
start_time = time.time()
tracemalloc.start()

# Read sequences from .fa files
collection_sequence = read_fasta('tester.sixteenth.fa')
reference_sequence = read_fasta('reference.fa')

if collection_sequence is None or reference_sequence is None:
    print("Error: One or both sequences could not be read.")
else:
    # Step 2: Preprocess reference
    reference_sequence = preprocess_reference(reference_sequence)
    
    # Step 3: Build the generalized suffix array
    gsa, isa, lcp, suffix_type = build_generalized_suffix_array(collection_sequence, reference_sequence)
    
    # Stop measuring time and space complexity
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Output the results
    '''print("Generalized Suffix Array (GSA):")
    for idx, stype in zip(gsa, suffix_type):
        print(f"{idx} ({stype})")
    
    print("\nInverted Suffix Array (ISA):", isa)
    print("Longest Common Prefix (LCP):", lcp)'''
    print(f"\nTime taken: {end_time - start_time:.5f} seconds")
    print(f"Current memory usage: {current / 1024 / 1024:.4f} MB")
    print(f"Peak memory usage: {peak / 1024 / 1024:.4f} MB")
