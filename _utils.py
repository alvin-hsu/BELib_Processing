from collections import namedtuple

REVCOMP_TT = str.maketrans('ATCGNWSMKRY',
                           'TAGCNWSKMYR')
Library = namedtuple('Library', ['df', 'lshmap'])
workers = set()


def revcomp(string):
    """ 
    Returns the reverse complement of the DNA sequence STRING.
    """
    return string.upper().translate(REVCOMP_TT)[::-1]


class KMerLSH(object):
    def __init__(self, entries, k=6):
        self.entries = []
        self.k = k 
        for idx, entry in enumerate(entries):
            s = set()
            for i in range(len(entry) - k): 
                s.add(entry[i:i+k])
            self.entries.append((idx, entry, s)) 

    def lookup(self, sequence, top_n=5, cutoff=30):
        kmers = []
        for i in range(len(sequence) - self.k):
            kmers.append(sequence[i:i+self.k])

        candidates = []
        for entry in self.entries:
            mismatches = 0 
            for kmer in kmers:
                if kmer not in entry[2]:
                    mismatches += 1
                    if mismatches > cutoff:
                        break
            else:
                candidates.append((entry, mismatches))

        candidates.sort(key=lambda x: x[1])
        return candidates[:top_n]


def poswise_match(query, reference, start=0, end=-1, stop=0):
    """ 
    Finds the REFERENCE string in the QUERY by sliding along QUERY starting from
    START and ending at END. If an alignment has a score <= STOP, then stops and
    returns immediately. Otherwise, returns the best index and its score.
    """
    assert len(query) >= len(reference)

    # end_idx is the last starting index to try
    if end < 0:
        end = len(query) + end
        end_idx = end - len(reference) + 1
    else:
        end_idx = end + 1

    best_score = len(reference)
    best_idx = None
    # We use memoryviews so Python doesn't make a new str every time we iterate
    query_b = memoryview(query.encode())
    ref_b = memoryview(reference.encode())
    for i in range(start, end_idx + 1):
        score = sum([c0 != c1 for c0, c1 in zip(ref_b, query_b[i:i + len(reference)])])
        if score < best_score:
            best_score = score
            best_idx = i
            if score <= stop:
                return i, score

    return best_idx, best_score


def make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = 0
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0
    return aligner

