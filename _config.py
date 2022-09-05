##########################
###  GENERAL SETTINGS  ###
##########################
# TODO: Set these to match your machine
N_WORKERS = 14  # Number of workers (Go with #cores - 2)

############################
###  ALIGNMENT SETTINGS  ###
############################
# Settings for the target read
TARGET_START = 0
TARGET_END = 56
TARGET_REVCOMP = True
SPACER_OFFSET = 21
# Settings for the gRNA read
GRNA_START = 18
GRNA_END = 42
GRNA_REVCOMP = False
# Quality cutoff for filtering reads
QUAL_CUTOFF = 28

#############################
###  PROCESSING SETTINGS  ###
#############################

