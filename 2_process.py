import logging
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy.stats import binom, f_oneway

from _config import *


def summarize(genotype_path: Path, lib_design: pd.DataFrame, name: str,
              error_pos_edits: Optional[np.ndarray] = None):
    """
    Generates a summary CSV file of the genotypes observed at GENOTYPE_PATH and saves it as
    `NAME_raw_stats.csv` if ERROR_POS_EDITS is None. Otherwise, will discard edits at each edit
    specified in ERROR_POS_EDITS and save the result as `NAME_corr_stats.csv`.
    """
    raw = (error_pos_edits is None)
    error_pos_edits = (np.zeros((len(lib_design), WINDOW_END - WINDOW_START + 1, len(BASE_EDITS)),
                                dtype=np.bool)
                       if error_pos_edits is None else error_pos_edits)
    stats = {'total_reads': np.zeros(len(lib_design)),
             'indel_reads': np.zeros(len(lib_design)),
             'unedited_reads': np.zeros(len(lib_design)),
             'baseedit_reads': np.zeros(len(lib_design)),
             'cbe_reads': np.zeros(len(lib_design)),
             'abe_reads': np.zeros(len(lib_design)),
             'poswise_be': {k: np.zeros((len(lib_design), WINDOW_END - WINDOW_START + 1))
                            for k in BASE_EDITS}
             }
    for lib_idx, row in lib_design.iterrows():
        total_reads = 0
        indel_reads = 0
        unedited_reads = 0
        baseedit_reads = 0
        unedited_target = row['original'][PS_OFFSET+WINDOW_START : PS_OFFSET+WINDOW_END+1]
        fixed_genotypes = []
        # Set the positions that can't possibly be base edited to nan
        for i, x in enumerate(unedited_target):
            for u, e in BASE_EDITS:
                if x != u:
                    stats['poswise_be'][u + e][lib_idx, i] = np.nan
        # Iterate over genotypes for each library member
        with (genotype_path / f'{lib_idx}.txt').open('r') as f:
            f.readline()  # First line of table is headers
            for line in f:
                genotype, count = line.strip().split(',')
                count = int(count)
                total_reads += count
                if genotype == 'INDEL':
                    indel_reads += count
                    continue
                genotype = genotype[PS_OFFSET + WINDOW_START : PS_OFFSET+WINDOW_END+1]
                # Fix genotype for random sequencing errors, assuming that "impossible" edits or
                # ambiguous positions are simply the unedited sequence.
                fixed_genotype = list(genotype)
                for i, (u, g) in enumerate(zip(unedited_target, genotype)):
                    if u != g and (u + g not in stats['poswise_be'] or
                                   error_pos_edits[lib_idx, i, BASE_EDITS.index(u + g)]):
                        fixed_genotype[i] = u
                fixed_genotype = ''.join(fixed_genotype)
                fixed_genotypes.append((fixed_genotype, count))
                # Record the stats for the fixed genotype
                if fixed_genotype == unedited_target:
                    unedited_reads += count
                else:
                    baseedit_reads += count
                has_cbe = False
                has_abe = False
                for i, (u, g) in enumerate(zip(unedited_target, fixed_genotype)):
                    if u != g:
                        edit = u + g
                        stats['poswise_be'][edit][lib_idx, i] += count
                        if edit == 'CT':
                            has_cbe = True
                        elif edit == 'AG':
                            has_abe = True
                if has_cbe:
                    stats['cbe_reads'][lib_idx] += count
                if has_abe:
                    stats['abe_reads'][lib_idx] += count
        # Compute average CA disequilibrium score and conditional probability
        c_pos = [i for i, x in enumerate(unedited_target, WINDOW_START)
                 if x == 'C' and BE_WIN_START <= i <= BE_WIN_END]
        a_pos = [i for i, x in enumerate(unedited_target, WINDOW_START)
                 if x == 'A' and BE_WIN_START <= i <= BE_WIN_END]
        ca_scores = []
        a_c_probs = []
        for c_i in c_pos:
            for a_i in a_pos:
                joint = 0
                c_total = 0
                a_total = 0
                for genotype, count in fixed_genotypes:
                    if genotype[c_i - WINDOW_START] == 'T' and genotype[a_i - WINDOW_START] == 'G':
                        joint += count
                    if genotype[c_i - WINDOW_START] == 'T':
                        c_total += count
                    if genotype[a_i - WINDOW_START] == 'G':
                        a_total += count
                if joint > 0 and a_total > 0 and c_total > 0:
                    ca_scores.append((abs(c_i - a_i), total_reads*joint / (a_total*c_total)))
                    a_c_probs.append((abs(c_i - a_i), joint / c_total))
        average_ca_scores = np.array([np.mean([x for d, x in ca_scores if d == dist])
                                      for dist in range(1, BE_WIN_END - BE_WIN_START+1)])
        average_a_c_probs = np.array([np.mean([x for d, x in a_c_probs if d == dist])
                                      for dist in range(1, BE_WIN_END - BE_WIN_START+1)])
        # Write statistics to report
        stats['total_reads'][lib_idx] = total_reads
        stats['indel_reads'][lib_idx] = indel_reads
        stats['unedited_reads'][lib_idx] = unedited_reads
        stats['baseedit_reads'][lib_idx] = baseedit_reads
        stats['ca_disequilibrium'][lib_idx, :] = average_ca_scores
        stats['a_given_c_prob'][lib_idx, :] = average_a_c_probs
        
    stats_df = {'total_reads': stats['total_reads'],
                'indel_reads': stats['indel_reads'],
                'unedited_reads': stats['unedited_reads'],
                'baseedit_reads': stats['baseedit_reads'],
                'cbe_reads': stats['cbe_reads'],
                'abe_reads': stats['abe_reads']}

    for i, pos in enumerate(range(WINDOW_START, WINDOW_END+1)):
        for u, e in BASE_EDITS:
            stats_df[f'pos_{pos}_{u}_to_{e}_count'] = stats['poswise_be'][u + e][:, i]
            stats_df[f'pos_{pos}_{u}_to_{e}_frac'] = stats_df[f'pos_{pos}_{u}_to_{e}_count'] / stats_df['total_reads']
    for dist in range(1, BE_WIN_END - BE_WIN_START+1):
        stats_df[f'ca_disequil_{dist}'] = stats['ca_disequilibrium'][:, dist-1]
        stats_df[f'a_given_c_prob_{dist}'] = stats['a_given_c_prob'][:, dist-1]

    stats_df = pd.DataFrame(stats_df)
    if raw:
        stats_df.to_csv(genotype_path / f'{name}_raw_stats.csv')
    else:
        stats_df.to_csv(genotype_path / f'{name}_corr_stats.csv')


def illumina_error_p(mut_cts, tot_cts, error_p=1e-3):
    """
    Returns the probability that a mutation present in each replicate with MUT_CTS reads out of
    TOT_CTS reads is due to Illumina error, rather than a real signal.
    """
    p = 1.0
    for mut_ct, tot_ct in zip(mut_cts, tot_cts):
        p = p * binom.sf(mut_ct, tot_ct, error_p)
    return p


def batch_correct(raw_summary_files: List[Path], lib_design: pd.DataFrame) -> np.ndarray:
    """
    Analyzes a batch of summary files for whether any batch effects or Illumina errors were
    observed. Returns an np.
    """
    summary_dfs = [pd.read_csv(x, index_col=0) for x in raw_summary_files]
    lib_size = len(lib_design)
    col_sets = [[f'pos_{pos}_{u}_to_{e}_frac' for u, e in BASE_EDITS] for pos in
                range(WINDOW_START, WINDOW_END + 1)]
    freqs = [np.array([[df[col_set].values.T for col_set in col_sets] for df in summary_dfs])]
    # Compute likelihood of batch effects, which result from uneven expansion of the cell library.
    p_batch_effect = np.zeros((lib_size, WINDOW_END - WINDOW_START + 1, len(BASE_EDITS)))
    for i in range(lib_size):
        for j, pos in enumerate(range(WINDOW_START, WINDOW_END + 1)):
            for k, _ in enumerate(BASE_EDITS):
                # Compute the probability that the editing results between these batches come from
                # the same distribution (one-way ANOVA). A low p-value indicates the presence of a
                # batch effect, which we should ignore.
                p_batch_effect[i, j, k] = f_oneway(*[a[:, j, k, i] for a in freqs])
    p_batch_flat = p_batch_effect.flatten()
    not_nan_inf = p_batch_flat[~(np.isnan(p_batch_flat) | np.isinf(p_batch_flat))]
    p_batch_cutoff = BATCH_FDR / not_nan_inf.size
    logging.info(f'Bonferonni-corrected p-cutoff for batch effects: {p_batch_cutoff}')
    logging.info(f'Minimum batch-effect p-value found: {not_nan_inf.min()}')
    logging.info(f'Number of batch effects found: {(not_nan_inf < p_batch_cutoff).sum()}')
    # Compute the likelihood that a base edit is simply due to Illumina sequencing error, rather
    # than due to an actual edit being present.
    p_illumina = np.ones((lib_size, WINDOW_END - WINDOW_START + 1, len(BASE_EDITS)))
    for i in range(lib_size):
        unedited = lib_design['original'][i][PS_OFFSET + WINDOW_START: PS_OFFSET + WINDOW_END + 1]
        total_reads = [df['total_reads'][i] for df in summary_dfs]
        for j, pos in enumerate(range(WINDOW_START, WINDOW_END + 1)):
            original_base = unedited[j]
            for k, (u, e) in enumerate(BASE_EDITS):
                if u == original_base:
                    mutation_cts = [df[f'pos_{pos}_{u}_to_{e}_count'][i] for df in summary_dfs]
                    p_illumina[i, j, k] = illumina_error_p(mutation_cts, total_reads)
    logging.info(f'Keeping {(p_illumina <= ILLUMINA_FDR).sum()} mutations')
    return (p_batch_effect < p_batch_cutoff) | (p_illumina > ILLUMINA_FDR)


def main(data_path: Path, lib_path: Path, exp_path: Path):
    lib_design = pd.read_csv(lib_path)
    exp_design = pd.read_csv(exp_path)
    exp_design['full_name'] = (exp_design['editor_name'] + '_' +
                               exp_design['cell_type'] + '_R' +
                               exp_design['rep'])
    # Generate raw summaries of each individual sequencing experiment
    for _, row in exp_design.iterrows():
        summarize(data_path / row['name'], lib_design, row['full_name'])
    # Batch-correct each group of experiments
    for editor in exp_design['editor_name'].unique():
        for cell_type in exp_design['cell_type'].unique():
            filtered = exp_design[(exp_design['editor_name'] == editor) &
                                  (exp_design['cell_type'] == cell_type)]
            if len(filtered) == 0:
                continue
            raw_summary_files = [data_path / row['name'] / f"{row['full_name']}_raw_stats.csv"
                                 for _, row in filtered.iterrows()]
            errors = batch_correct(raw_summary_files)
            for _, row in filtered.iterrows():
                summarize(data_path / row['name'], lib_design, row['full_name'], errors)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('data_dir', type=str,
                        help='Directory with *_genotypes/ folders')
    parser.add_argument('lib_design', type=str,
                        help='Library design CSV file (see ./lib_designs/README.md)')
    parser.add_argument('exp_design', type=str,
                        help='Experimental design CSV file (For example, see ./MN_exp_design.csv)')
    args = parser.parse_args()
    data_path = Path(args.data_dir).resolve()
    lib_path = Path(args.lib_design).resolve()
    exp_path = Path(args.exp_design).resolve()
    assert data_path.is_dir(), 'Missing data path'
    assert lib_path.is_file(), 'Missing library desisgn CSV'
    assert exp_path.is_file(), 'Missing experimental desisgn CSV'
    main(data_path, lib_path, exp_path)
