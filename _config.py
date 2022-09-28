##########################
###  GENERAL SETTINGS  ###
##########################
# Set these to match your machine!
N_WORKERS = 14  # Number of workers (Go with #cores - 2)

############################
###  ALIGNMENT SETTINGS  ###
############################
# Settings for the target read
TARGET_START = 0
TARGET_LEN = 56
TARGET_REVCOMP = True
# Settings for the gRNA read
GRNA_START = 18
GRNA_END = 42
GRNA_REVCOMP = False
# Quality cutoff for filtering reads
QUAL_CUTOFF = 28

#############################
###  PROCESSING SETTINGS  ###
#############################
PS_OFFSET = 21
BASE_EDITS = ['AG', 'CA', 'CG', 'CT']
WINDOW_START = -9
WINDOW_END = 20
BE_WIN_START = 3
BE_WIN_END = 8
MIN_READCOUNT = 1000
BATCH_FDR = 0.005
ILLUMINA_FDR = 0.05

###########################
###  ANALYSIS SETTINGS  ###
###########################
ANALYSIS_NAMES = ['satmut_6mer', 'satmut_4mer']
POS_COUNT_COLS = [f'pos_{pos}_{u}_to_{e}_count'
                  for pos in range(WINDOW_START, WINDOW_END+1)
                  for u, e in BASE_EDITS]
SUM_COLS = (['total_reads', 'indel_reads', 'baseedit_reads', 'cbe_reads', 'abe_reads', 'crbe_reads',
             'not_cbe_reads', 'pure_reads', 'unedited_reads'] + POS_COUNT_COLS)
LOGIT_EPS = 0.001
LOGO_START = 4
LOGO_END = 8
LOGO_LSIZE = 3
LOGO_RSIZE = 3
