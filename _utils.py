from collections import namedtuple
from datetime import datetime
import logging
from pathlib import Path
from subprocess import check_output
from typing import Tuple

from Bio.Align import PairwiseAligner


Library = namedtuple('Library', ['df', 'lshmap'])


REVCOMP_TT = str.maketrans('ATCGNWSMKRY',
                           'TAGCNWSKMYR')
def revcomp(string):
    """ 
    Returns the reverse complement of the DNA sequence STRING.
    """
    return string.upper().translate(REVCOMP_TT)[::-1]


class KMerLSH(object):
    """
    A locality-sensitive hashing table.
    """
    def __init__(self, entries, k=6):
        self.entries = []
        self.k = k 
        for idx, entry in enumerate(entries):
            s = set()
            for i in range(len(entry) - k): 
                s.add(entry[i:i+k])
            self.entries.append((idx, entry, s)) 

    def lookup(self, sequence, top_n=5, cutoff=30) -> Tuple[str, int]:
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


def dt_str() -> str:
    """
    Returns a string with the current datetime.
    """
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def git_commit() -> str:
    """
    Returns a string with the first digits of the commit hash.
    """
    git_command = ['git', 'rev-parse', '--short', 'HEAD']
    return check_output(git_command).decode('ascii').strip()


def configure_logger(log_path: Path):
    """
    Configures the Python `logging` package to print pretty logging messages that will be stored at 
    LOG_PATH/alignment_log_YYMMDD_HHMMSS.txt
    """
    log_path = log_path / f'alignment_log_{dt_str()}.txt'
    log_path.touch()
    logging.getLogger().setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[logging.FileHandler(log_path),
                  logging.StreamHandler()]
    )
    logging.info(f'Saving logs to {str(log_path)}')
    logging.info(f'Running analysis with commit ID {git_commit()}')


def poswise_match(query: str, reference: str, start=0, end=-1, stop=0) -> Tuple[int, int]:
    """ 
    Finds the REFERENCE string in the QUERY by sliding along QUERY starting from START and ending at
    END. If an alignment has a score <= STOP, then stops and returns immediately. Otherwise, returns
    the best index and its score.
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


def make_aligner() -> PairwiseAligner:
    """
    Returns a PairwiseAligner object with the parameters that we typically use for analysis.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = 0
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0
    return aligner

