def dt_str() -> str:
    """
    Returns a string with the current datetime.
    """
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def git_commit() -> str:
    """
    Returns a string with the first digits of the commit hash.
    """
    from subprocess import check_output
    git_command = ['git', 'rev-parse', '--short', 'HEAD']
    return check_output(git_command).decode('ascii').strip()


def configure_logger(log_path: str):
    """
    Configures the Python `logging` package to print pretty logging messages that will be stored at 
    LOG_PATH/[script_name]_log_YYMMDD_HHMMSS.txt
    """
    import logging
    from pathlib import Path
    from sys import argv
    script_name = Path(argv[1]).stem
    log_path = Path(log_path).resolve() / f'{script_name}_log_{dt_str()}.txt'
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


def configure_mpl():
    """
    Configures matplotlib to produce PDFs compatible with Illustrator.
    """
    import os
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    fonts = fm.findSystemFonts(fontpaths=[os.path.join(plt.rcParams['datapath'], 'fonts/ttf')])
    for font in fonts:
        fm.fontManager.addfont(font)
    matplotlib.use('pgf')
    plt.rcParams['ps.useafm'] = True
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rcParams['savefig.transparent'] = True
    plt.rcParams['pgf.texsystem'] = 'xelatex'


REVCOMP_TT = str.maketrans('ATCGNWSMKRY',
                           'TAGCNWSKMYR')
def revcomp(string):
    """ 
    Returns the reverse complement of the DNA sequence STRING.
    """
    return string.upper().translate(REVCOMP_TT)[::-1]


def make_aligner():
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = 0
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0
    return aligner


BASES = list('ACGT')
