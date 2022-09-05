from collections import defaultdict
import gzip
import logging
from multiprocessing import Lock, Manager, Process, Pool, Queue, Value
import os
import pickle as pkl
from sys import argv
from time import sleep

import numpy as np
import pandas as pd
from pexpect import spawn
from Bio.Align import PairwiseAligner

from _config import *
from _utils import *

LIB_PKL = 'library_12kChar_AH_fmt.pkl'

def process_reads(q, locks, genotypes, total_q_pass, total_workers_done):

    with open(LIB_PKL, 'rb') as f:
        lib_df, lib_lsh = pkl.load(f)

    q_pass = 0
    while True:
        read_pair = q.get()
        # When all reads have been pushed, master process puts 'stop' in queue
        # for each worker to signal that there are no more reads.
        if read_pair == 'stop':
            break
        try:
            read1, qual1, read2, qual2 = read_pair
        except ValueError:
            continue
        # Check the quality of the relevant regions to filter low-quality reads
        tgt_qual = qual1[TARGET_START:TARGET_END]
        grna_qual = qual2[GRNA_START:GRNA_END]
        qual_bytes = (tgt_qual + grna_qual).encode()
        mean_q = np.mean(np.frombuffer(qual_bytes, dtype=np.uint8).astype(np.float32)) - 33
        if mean_q < QUAL_CUTOFF:
            continue
        q_pass += 1
        target_read = read1[TARGET_START:TARGET_END]
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


def write_genotypes(out_dir, idx, d):
    with open(out_dir + f'/{idx}.txt', 'w+') as f:
        f.write(',count\n')
        for k, v in d.items():
            f.write(f'{k},{v}\n')


def main(prefix):
    r1_fname = prefix + '_R1_001.fastq.gz'
    r2_fname = prefix + '_R2_001.fastq.gz'
    out_dir = prefix + '_out'
    os.makedirs(out_dir, exist_ok=True)
    # TODO: Change this to csv and make the pkl if it doesn't exist
    with open(LIB_PKL, 'rb') as f:
        lib_df, lib_lsh = pkl.load(f)
    # Synchronization variables
    q = Queue()  # Queue of reads to process
    total_q_pass = Value('L', 0)  # Number of reads passing filter
    total_workers_done = Value('I', 0)  # Number of workers that have finished
    with Manager() as manager:
        locks = [Lock() for _ in range(len(lib_df))]  # A lock to prevent race conditions
        genotypes = [manager.dict() for _ in range(len(lib_df))]  # Genotypes for each library member
        workers = [Process(target=process_reads,
                           args=(q, locks, genotypes, total_q_pass, total_workers_done))
                   for _ in range(N_WORKERS)]
        print(f'Spawning {N_WORKERS} workers...', end='')
        for proc in workers:
            proc.start()
            sleep(0.1)
        print(f'done')
        with gzip.open(r1_fname, 'r') as f1, gzip.open(r2_fname, 'r') as f2:
            for i, (line1, line2) in enumerate(zip(f1, f2)):
                if i % 4 == 1:
                    read1 = line1.decode().strip()
                    read2 = line2.decode().strip()
                elif i % 4 == 3:
                    qual1 = line1.decode().strip()
                    qual2 = line2.decode().strip()
                    q.put(((read1, qual1, read2, qual2)))
                if i % 4000000 == 0 and i > 0:
                    print(f'Added {i // 4} reads to the queue. Remaining: {q.qsize()}')
        for _ in range(N_WORKERS):
            q.put('stop')

        # Wait for workers to finish
        last_remaining = q.qsize()
        print(f'Added {i // 4} reads to the queue. Remaining: {last_remaining} reads')
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
            print(f'Remaining: {now_remaining} reads [{reads_per_min} reads/min]')
            last_remaining = now_remaining
        print('All workers done aligning!')
        # Stop all of the worker processes
        for proc in workers:
            proc.join()
        # Make new worker processes that process each genotype dict.
        with Pool(N_WORKERS) as pool:
            pool.starmap(write_genotypes, [(out_dir, i, genotypes[i]) for i in range(len(lib_df))])


if __name__ == '__main__':
    # TODO: Make this argparse
    from sys import argv
    prefix = argv[1]
    main(prefix)

