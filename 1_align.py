import gzip
from itertools import zip_longest
import logging
from multiprocessing import Lock, Manager, Process, Pool, Queue, Value
from pathlib import Path
import pickle as pkl
from time import sleep
from typing import List, Dict

import numpy as np
import pandas as pd
from pexpect import spawn

from _config import *
from _utils import *


def process_reads(queue: Queue, pkl_fname: str,
                  locks: List[Lock], genotypes: List[Dict[str, int]],
                  total_q_pass: int, total_workers_done: int):
    """
    Code that a worker process should follow during alignment. Note that LOCKS,
    GENOTYPES, TOTAL_Q_PASS, and TOTAL_WORKERS_DONE are not actually the types
    that their typehints would imply, but rather thread-safe implementations
    that work with the same API as suggested by the typehints.
    """
    with open(pkl_fname, 'rb') as f:
        lib_df, lib_lsh = pkl.load(f)

    q_pass = 0
    while True:
        read_pair = queue.get()
        # When all reads have been pushed, master process puts 'stop' in queue
        # for each worker to signal that there are no more reads.
        if read_pair == 'stop':
            break
        try:
            read1, qual1, read2, qual2 = read_pair
        except ValueError:
            continue
        # Check the quality of the relevant regions to filter low-quality reads
        tgt_qual = qual1[TARGET_START:TARGET_START+TARGET_LEN]
        grna_qual = qual2[GRNA_START:GRNA_END]
        qual_bytes = (tgt_qual + grna_qual).encode()
        mean_q = np.mean(np.frombuffer(qual_bytes, dtype=np.uint8).astype(np.float32)) - 33
        if mean_q < QUAL_CUTOFF:
            continue
        q_pass += 1
        target_read = read1[TARGET_START:TARGET_START+TARGET_LEN]
        target_read = revcomp(target_read) if TARGET_REVCOMP else target_read
        grna_read = read2[GRNA_START:GRNA_END]

        # Find likely target matches with locality-sensitive hashing
        tgt_matches = lib_lsh.lookup(target_read, top_n=5)
        if not tgt_matches:
            continue
        grna_matches = []
        for candidate in tgt_matches:
            idx = candidate[0][0]
            pos, score = poswise_match(grna_read, lib_df.iloc[idx]['gRNA'],
                                       start=18, end=23, stop=2)
            if pos is not None:
                grna_matches.append((idx, score))
        if not grna_matches:
            continue
        grna_matches.sort(key=lambda x: x[1])
        best_idx = grna_matches[0][0]
        best_target = lib_df.iloc[best_idx]['original']

        # Align the best candidate with the read to check for indels
        alignment = aligner.align(best_target, target_read)[0]
        tgt_align, align_str, read_align, _ = str(alignment).split('\n')
        if '-' in align_str:
            genotype = 'INDEL'
        else:
            genotype = read_align

        locks[best_idx].acquire()
        if genotype in genotypes[best_idx]:
            genotypes[best_idx][genotype] += 1
        else:
            genotypes[best_idx][genotype] = 1
        locks[best_idx].release()

    # Reports number of passed reads
    with total_q_pass.get_lock():
        total_q_pass.value += q_pass
    # Signals to master process that worker process is finished
    with total_workers_done.get_lock():
        total_workers_done.value += 1


def write_genotypes(out_path: Path, idx: int, d: Dict[str, int]):
    """
    Write genotypes from D to OUT_PATH/IDX.txt.
    """
    with (out_path / f'{idx}.txt').open('w+') as f:
        f.write('genotype,count\n')
        for k, v in d.items():
            f.write(f'{k},{v}\n')


def main(tgt_path: Path, grna_path: Path, lib_path: Path):
    """
    Pushes paired-end reads from TGT_PATH and GRNA_PATH to a multiprocessing queue, where workers
    (also spawned by this function) will align the results. Then maps the library members over
    new workers for processing the genotypes assigned to each library member.
    """
    # Make an output folder, if it doesn't already exist.
    for i, (t, g) in enumerate(zip_longest(str(tgt_path), str(grna_path))):
        if t != g:
            break
    prefix = str(tgt_path)[:i].rsplit('_', 1)[0]
    out_path = Path(prefix + '_genotypes').resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    configure_logger(out_path)
    # Check to see if the library information is already cached as a pickle. If not, create it.
    pkl_path = lib_path.with_suffix('.pkl')
    if pkl_path.is_file():
        with pkl_path.open('rb') as f:
            lib_df, lib_lsh = pkl.load(f)
    else:
        lib_df = pd.read_csv(lib_path)
        library = Library(lib_df, KMerLSH(lib_df['original']))
        with pkl_path.open('wb+') as f:
            pkl.dump(library, f)
        lib_df, lib_lsh = library

    # Set up synchronization primitives
    q = Queue()                         # Queue of reads to process
    total_q_pass = Value('L', 0)        # Number of reads passing filter
    total_workers_done = Value('I', 0)  # Number of workers that have finished
    with Manager() as manager:
        locks = [Lock() for _ in range(len(lib_df))]              # Locks to prevent race conditions
        genotypes = [manager.dict() for _ in range(len(lib_df))]  # Genotypes for each member
        
        logging.info(f'Spawning {N_WORKERS} workers...')
        workers = [Process(target=process_reads,
                           args=(q, locks, genotypes, total_q_pass, total_workers_done))
                   for _ in range(N_WORKERS)]
        for proc in workers:
            proc.start()
            sleep(0.1)
        logging.info(f'Finished')
        with gzip.open(str(tgt_path), 'r') as f1, gzip.open(str(grna_path), 'r') as f2:
            for i, (line1, line2) in enumerate(zip(f1, f2)):
                if i % 4 == 1:
                    read1 = line1.decode().strip()
                    read2 = line2.decode().strip()
                elif i % 4 == 3:
                    qual1 = line1.decode().strip()
                    qual2 = line2.decode().strip()
                    q.put(((read1, qual1, read2, qual2)))
                if i % 4000000 == 0 and i > 0:
                    logging.info(f'Added {i // 4} reads to the queue. Remaining: {q.qsize()}')
        for _ in range(N_WORKERS):
            q.put('stop')

        # Wait for workers to finish
        last_remaining = q.qsize()
        logging.info(f'Added {i // 4} reads to the queue. Remaining: {last_remaining} reads')
        cons_0 = 0  # Consecutive minutes the number of reads hasn't decreased. Sometimes one
                    # worker will mysteriously fail to stop and deadlock the program. This
                    # allows the program to finish running if that happens.
        while total_workers_done.value < N_WORKERS:
            sleep(60)
            now_remaining = q.qsize()
            reads_per_min = last_remaining - now_remaining
            if reads_per_min == 0:
                cons_0 += 1
                if cons_0 == 10:
                    break
            else:
                cons_0 = 0
            logging.info(f'Remaining: {now_remaining} reads [{reads_per_min} reads/min]')
            last_remaining = now_remaining
        logging.info('All workers done aligning!')
        # Stop all of the worker processes
        for proc in workers:
            proc.join()
        # Make new worker processes that process each genotype dict.
        with Pool(N_WORKERS) as pool:
            pool.starmap(write_genotypes, [(out_path, i, genotypes[i]) for i in range(len(lib_df))])


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('tgt_fname', type=str,
                        help='Read .fastq.gz file containing the target (usually R1)')
    parser.add_argument('grna_fname', type=str,
                        help='Read .fastq.gz file containing the gRNA (usually R2)')
    parser.add_argument('lib_design', type=str,
                        help='Library design CSV file (see ./lib_designs/README.md)')
    args = parser.parse_args()
    tgt_path = Path(args.tgt_fname).resolve()
    grna_path = Path(args.grna_fname).resolve()
    lib_path = Path(args.lib_design).resolve()
    assert tgt_path.is_file(), 'Missing target .fastq.gz (usually is R1)'
    assert grna_path.is_file(), 'Missing gRNA .fastq.gz (usually is R2)'
    assert lib_path.is_file(), 'Missing library desisgn CSV'
    main(tgt_path, grna_path, lib_path)

