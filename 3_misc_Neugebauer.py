import logging
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy.odr as odr

from _config import *
from _utils import *


PAIRS = [(('TadCBEa', 'mES'), ('TadCBEa_V106W', 'mES')),
         (('TadCBEb', 'mES'), ('TadCBEb_V106W', 'mES')),
         (('TadCBEc', 'mES'), ('TadCBEc_V106W', 'mES')),
         (('TadCBEd', 'mES'), ('TadCBEd_V106W', 'mES'))]


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


def scatter_with_tls(xs: np.ndarray, ys: np.ndarray, pdf: PdfPages,
                     title: str = '', xlabel: str = '', ylabel: str = ''):
    lr = odr.Model(lambda b, x: b[0]*x)
    data = odr.Data(xs, ys)
    tls = odr.ODR(data, lr, beta0=[1.])
    output = tls.run()
    output.pprint()
    b1 = output.beta[0]
    plt.figure(figsize=(4, 4), dpi=50)
    plt.scatter(xs, ys, c='b', s=3, alpha=0.1)
    plt.plot([0, 1], [0, 1], c='r', label='y = x')
    plt.plot([0, 1], [0, b1], '--', c='r', label=f'TLS line\n(y = {b1:.2f}x)')
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    pdf.savefig(bbox_inches='tight')
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
    v106w_pdf = PdfPages(out_path / 'v106w_plots.pdf')
    for (editor_1, cell_type_1), (editor_2, cell_type_2) in PAIRS:
        name_1 = f'{editor_1}_{cell_type_1}'
        logging.info(f'Processing {name_1}')
        filtered_1 = exp_design[(exp_design['editor_name'] == editor_1) &
                                (exp_design['cell_type'] == cell_type_1)]
        summary_files_1 = [(data_path /
                            f"{row['name']}_genotypes" /
                            f"{row['full_name']}_corr_stats.csv")
                           for _, row in filtered_1.iterrows()]
        df1, _ = combine_stats(summary_files_1)
        df1 = df1.join(lib_design)
        df1 = df1[df1['name'].isin(ANALYSIS_NAMES) & (df1['total_reads'] >= MIN_READCOUNT)]

        name_2 = f'{editor_2}_{cell_type_2}'
        logging.info(f'Processing {name_2}')
        filtered_2 = exp_design[(exp_design['editor_name'] == editor_2) &
                                (exp_design['cell_type'] == cell_type_2)]
        summary_files_2 = [(data_path /
                            f"{row['name']}_genotypes" /
                            f"{row['full_name']}_corr_stats.csv")
                           for _, row in filtered_2.iterrows()]
        df2, _ = combine_stats(summary_files_2)
        df2 = df2.join(lib_design)
        df2 = df2[df2['name'].isin(ANALYSIS_NAMES) & (df2['total_reads'] >= MIN_READCOUNT)]

        joined = df1.join(df2, how='inner', lsuffix='_1', rsuffix='_2').dropna()
        scatter_with_tls(joined['pos_6_C_to_T_frac_1'].values, joined['pos_6_C_to_T_frac_2'].values,
                         v106w_pdf, 'C-to-T editing efficiency at protospacer position 6',
                         f'{name_1}_efficiency', f'{name_2}_efficiency')


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
