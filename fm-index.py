import time
import tracemalloc

# Provided functions and classes
def suffixArray(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    return list(map(lambda x: x[1], satups)) # extract, return just offsets

def bwtFromSa(t, sa=None):
    bw = []
    dollarRow = None
    if sa is None:
        sa = suffixArray(t)
    for si in sa:
        if si == 0:
            dollarRow = len(bw)
            bw.append('$')
        else:
            bw.append(t[si-1])
    return ''.join(bw), dollarRow 

class FmCheckpoints(object):
    
    def __init__(self, bw, cpIval=4):
        self.cps = {}        # checkpoints
        self.cpIval = cpIval # spacing between checkpoints
        tally = {}           # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        for i, c in enumerate(bw):
            tally[c] += 1 
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])
    
    def rank(self, bw, c, row):
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc

class FmIndex():
    
    @staticmethod
    def downsampleSuffixArray(sa, n=4):
        ssa = {}
        for i, suf in enumerate(sa):
            if suf % n == 0:
                ssa[i] = suf
        return ssa
    
    def __init__(self, t, cpIval=4, ssaIval=4):
        if t[-1] != '$':
            t += '$' # add dollar if not there already
        # Get BWT string and offset of $ within it
        sa = suffixArray(t)
        self.bwt, self.dollarRow = bwtFromSa(t, sa)
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.slen = len(self.bwt)
        self.cps = FmCheckpoints(self.bwt, cpIval)
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count    
    
    def occurrences(self, p):
        left, right = 0, self.slen - 1
        for i in range(len(p)-1, -1, -1):
            left = self.first[p[i]] + self.cps.rank(self.bwt, p[i], left-1)
            right = self.first[p[i]] + self.cps.rank(self.bwt, p[i], right) - 1
            if right < left:
                return []

        occurrence_range = (left, right)
        if occurrence_range is None:
            return []

        left, right = occurrence_range
        offsets = []
        for i in range(left, right+1):
            row = i
            steps = 0
            while True:
                c = self.bwt[row]
                row = self.first[c] + self.cps.rank(self.bwt, c, row) - 1
                steps += 1
                if row in self.ssa:
                    offsets.append(self.ssa[row] + steps)
                    break
        return offsets

# Additional functions for reading FASTA files and measuring time/space complexity
def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def build_generalized_suffix_array(collection, reference):
    combined_string = collection + '$' + reference + '#'
    sa = suffixArray(combined_string)
    # Determine which suffixes belong to the collection or the reference
    suffix_type = ['C' if i < len(collection) else 'R' for i in sa]
    return sa, suffix_type

# Measure time and space complexity
start_time = time.time()
tracemalloc.start()

# Read sequences from .fa files
collection_sequence = read_fasta('tester.sixteenth.fa')
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
