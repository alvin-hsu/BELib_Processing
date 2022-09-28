import logging
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import gmean, pearsonr
from sklearn.linear_model import Ridge
import logomaker

from _config import *
from _utils import *


def combine_stats(files: List[Path]) -> Tuple[pd.DataFrame, List[pd.DataFrame]]:
    stats_dfs = [pd.read_csv(p) for p in files]
    combined_df = {}
    for col in SUM_COLS:
        combined_df[col] = stats_dfs[0]
        for df in stats_dfs[1:]:
            combined_df[col] = combined_df[col] + df[col]
        for pos in range(WINDOW_START, WINDOW_END + 1):
            for u, e in BASE_EDITS:
                combined_df[f'pos_{pos}_{u}_to_{e}_frac'] = \
                    combined_df[f'pos_{pos}_{u}_to_{e}_count'] / combined_df['total_reads']
    return pd.DataFrame(combined_df), stats_dfs


def plot_rep_corrs(df1: pd.DataFrame, df2: pd.DataFrame, pdf: PdfPages,
                   edit: str, name: str, r1: int, r2: int):
    cols = [f'pos_{pos}_{edit[0]}_to_{edit[1]}_frac' for pos in range(WINDOW_START, WINDOW_END + 1)]
    df = df1[cols].join(df2[cols], how='inner', lsuffix=f'_{r1}', rsuffix=f'_{r2}')
    poswise_reps = pd.DataFrame()
    for pos in range(WINDOW_START, WINDOW_END + 1):
        poswise_reps = pd.concat((poswise_reps, pd.DataFrame(
            {f'rep_{r1}': df[f'pos_{pos}_{edit[0]}_to_{edit[1]}_frac_{r1}'],
             f'rep_{r2}': df[f'pos_{pos}_{edit[0]}_to_{edit[1]}_frac_{r2}']})))
    plt.figure(figsize=(4, 4), dpi=150)
    plt.scatter(poswise_reps['rep1'], poswise_reps['rep2'], s=1, alpha=0.1, c='b', rasterized=True)
    r = pearsonr(poswise_reps['rep1'], poswise_reps['rep2'])[0]
    plt.plot([0, 1], [0, 1], label=f'Pearson R = {r:.3f}', c='r')
    plt.xlabel(f'{edit[0]}-to-{edit[1]} activity (rep. {r1})')
    plt.ylabel(f'{edit[0]}-to-{edit[1]} activity (rep. {r2})')
    plt.title(f'Replicate consistency for {name}')
    plt.legend()
    pdf.savefig(bbox_inches='tight')
    plt.close()


def log_aggregate_stats(df: pd.DataFrame):
    logging.info(f"Average fraction of reads edited: "
                 f"{(df['baseedit_reads'] / df['total_reads']).mean()}")
    logging.info(f"Average fraction of reads with any C>T editing: "
                 f"{(df['cbe_reads'] / df['total_reads']).mean()}")
    logging.info(f"Average fraction of reads with any A>G editing: "
                 f"{(df['abe_reads'] / df['total_reads']).mean()}")
    ct_activity = [df[f'pos_{idx}_C_to_T_frac'].dropna().mean()
                   for idx in range(WINDOW_START, WINDOW_END + 1)]
    ag_activity = [df[f'pos_{idx}_A_to_G_frac'].dropna().mean()
                   for idx in range(WINDOW_START, WINDOW_END + 1)]
    max_ct = max(ct_activity)
    max_ag = max(ag_activity)
    c_to_a = [c / a for c, a in zip(ct_activity, ag_activity) if c >= 0.3 * max_ct]
    logging.info(f"Geometric mean of C>T:A>G ratio: {gmean(c_to_a)}")
    logging.info(f"Maximum average C>T editing activity: {max_ct}")
    logging.info(f"Maximum average A>G editing activity: {max_ag}")


def plot_ct_ag_windows(df: pd.DataFrame, pdf: PdfPages, name: str):
    ct_activity = [df[f'pos_{idx}_C_to_T_frac'].dropna().mean()
                   for idx in range(WINDOW_START, WINDOW_END + 1)]
    ag_activity = [df[f'pos_{idx}_A_to_G_frac'].dropna().mean()
                   for idx in range(WINDOW_START, WINDOW_END + 1)]
    max_ct = max(ct_activity)
    max_ag = max(ag_activity)
    in_ct_window = [(x / max_ct) >= 0.3 for x in ct_activity]
    in_ag_window = [(x / max_ag) >= 0.3 for x in ag_activity]
    # Plotting
    plt.figure(figsize=(5, 4), dpi=50)
    ax = plt.gca()
    # Plot A>G and C>T activity windows
    plt.plot(range(WINDOW_START, WINDOW_END + 1), ag_activity, 'o-', c='r', label='A>G editing')
    plt.plot(range(WINDOW_START, WINDOW_END + 1), ct_activity, 'o-', c='b', label='C>T editing')
    # Indicate the protospacer positions
    plt.axvline(0.5, c='k')
    plt.axvline(20.5, c='k')
    # Draw the window
    plt.axhline(0.3 * max_ct, linestyle='--', c='b', label='CBE window cutoff')
    for i, (c, a) in enumerate(zip(in_ct_window, in_ag_window), WINDOW_START):
        ax.add_patch(Rectangle((i - 0.5, 0), 1, c, facecolor='#ddddff80'))
    plt.xlim(-10, 20)
    plt.ylim(0, 1)
    plt.xlabel('Protospacer position')
    plt.ylabel('Average fraction of bases converted')
    plt.title(name)
    plt.legend()
    pdf.savefig()
    plt.close()


def plot_seq_logo(df: pd.DataFrame, pdf: PdfPages, name: str, edit: str):
    def stabilized_logit(x):
        return np.log((x + LOGIT_EPS) / (1 - x + LOGIT_EPS))

    # Scrape all contexts from library
    contexts = []
    ys = []
    for pos in range(LOGO_START, LOGO_END + 1):
        col_name = f'pos_{pos}_{edit[0]}_to_{edit[1]}_frac'
        filtered = df.loc[df[col_name].dropna().index]
        lstart, lend = PS_OFFSET + pos - LOGO_LSIZE, PS_OFFSET + pos
        rstart, rend = PS_OFFSET + pos + 1, PS_OFFSET + pos + LOGO_RSIZE + 1
        l_ctx = filtered['original'].str.slice(lstart, lend).tolist()
        r_ctx = filtered['original'].str.slice(rstart, rend).tolist()
        contexts = contexts + [left + right for left, right in zip(l_ctx, r_ctx)]
        ys = ys + filtered[col_name].apply(stabilized_logit).tolist()
    num_sites = len(contexts)
    # Make one-hot encoding matrices
    full_x = np.zeros((num_sites, 4 * (LOGO_LSIZE + LOGO_RSIZE)))
    full_y = np.array(ys)
    for i, context in enumerate(contexts):
        for j, base in enumerate(context):
            full_x[i, 4 * j + BASES.index(base)] = 1
    # Train/Test split
    full_idx = np.arange(num_sites, dtype=np.uint32)
    np.random.shuffle(full_idx)
    data_split = int(0.8 * num_sites)
    train_idx, test_idx = full_idx[:data_split], full_idx[data_split:]
    train_x, train_y = full_x[train_idx, :], full_y[train_idx]
    test_x, test_y = full_x[test_idx, :], full_y[test_idx]
    # Ridge regression on logits
    model = Ridge(1e-5, solver='svd')
    model.fit(train_x, train_y)
    test_r2 = model.score(test_x, test_y)
    coef = model.coef_.reshape((LOGO_LSIZE + LOGO_RSIZE), 4)
    target = np.array([[0, 0, 0, 0]], dtype=np.float32)
    target[0, BASES.index(edit[0])] = 1
    logo_data = np.concatenate([coef[:LOGO_LSIZE, :],
                                target,
                                coef[LOGO_LSIZE:, :]], axis=0)
    logo_df = pd.DataFrame(columns=['A', 'C', 'G', 'T'],
                           index=list(range(-LOGO_LSIZE, LOGO_RSIZE + 1)),
                           data=logo_data)
    color_scheme = {'A': [1.0, 0.2, 0.2],
                    'C': [0.2, 0.2, 1.0],
                    'G': [1.0, 0.8, 0.2],
                    'T': [0.2, 1.0, 0.2]}
    alpha = min(0.5, np.sqrt(test_r2)) / 0.5
    plt.figure(figsize=(0.75 * (LOGO_LSIZE + LOGO_RSIZE + 1), 4), dpi=50)
    ax = plt.gca()
    logo = logomaker.Logo(logo_df, ax=ax, flip_below=False, color_scheme=color_scheme,
                          font_name='Nimbus Sans', baseline_width=0.1, alpha=alpha)
    logo.glyph_df.loc[0, edit[0]].set_attributes(floor=-0.5, ceiling=0.5, color=[0.75, 0.75, 0.75])
    plt.ylabel('Logistic regression weights')
    plt.xlabel(f'Position relative to {edit[0]}')
    plt.xticks(list(range(-LOGO_LSIZE, 0)) + list(range(1, LOGO_RSIZE + 1)))
    plt.title(f'{name} (Pos. {LOGO_START}-{LOGO_END} [N={train_idx.size}]\n'
              f'Test R = {np.sqrt(test_r2):.4f})')
    logo.draw()
    pdf.savefig()
    plt.close()


def main(data_path: Path, lib_path: Path, exp_path: Path):
    configure_logger(str(data_path))
    configure_mpl()
    lib_design = pd.read_csv(lib_path)
    exp_design = pd.read_csv(exp_path)
    exp_design['full_name'] = (exp_design['editor_name'] + '_' +
                               exp_design['cell_type'] + '_R' +
                               exp_design['rep'])
    out_path = data_path / 'figures'
    rep_corr_pdf = PdfPages(out_path / 'replicate_correlation.pdf')
    window_pdf = PdfPages(out_path / 'windows.pdf')
    cbe_logo_pdf = PdfPages(out_path / 'cbe_logos.pdf')
    abe_logo_pdf = PdfPages(out_path / 'abe_logos.pdf')
    for editor in exp_design['editor_name'].unique():
        for cell_type in exp_design['cell_type'].unique():
            filtered = exp_design[(exp_design['editor_name'] == editor) &
                                  (exp_design['cell_type'] == cell_type)]
            if len(filtered) == 0:
                continue
            summary_files = [(data_path /
                              f"{row['name']}_genotypes" /
                              f"{row['full_name']}_corr_stats.csv")
                             for _, row in filtered.iterrows()]
            name = f'{editor}_{cell_type}'
            logging.info(f'Processing {name}')
            combined_df, ind_dfs = combine_stats(summary_files)
            ind_dfs = [x.join(lib_design) for x in ind_dfs]
            ind_dfs = [x[x['name'].isin(ANALYSIS_NAMES) & x['total_reads'] >= MIN_READCOUNT]
                       for x in ind_dfs]
            combined_df = combined_df.join(lib_design)
            combined_df = combined_df[combined_df['name'].isin(ANALYSIS_NAMES) &
                                      (combined_df['total_reads'] >= MIN_READCOUNT)]
            if len(summary_files) >= 2:
                for i in range(len(filtered)):
                    df_i = ind_dfs[i]
                    idx_i = filtered['rep'].iloc[i]
                    for j in range(i + 1, len(filtered)):
                        df_j = ind_dfs[j]
                        idx_j = filtered['rep'].iloc[j]
                        plot_rep_corrs(df_i, df_j, rep_corr_pdf, 'CT', name, idx_i, idx_j)
            log_aggregate_stats(combined_df, lib_design)
            plot_ct_ag_windows(combined_df, window_pdf, name)
            plot_seq_logo(combined_df, cbe_logo_pdf, name, 'CT')
            plot_seq_logo(combined_df, abe_logo_pdf, name, 'AG')
    rep_corr_pdf.close()
    window_pdf.close()
    cbe_logo_pdf.close()
    abe_logo_pdf.close()


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
    main()
