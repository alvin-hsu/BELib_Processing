# Base Editor Library Processing
Code for processing data from base editing screens performed with the libraries described in Arbab, Shen, et al. 2020. The code in this repository was used to process data for Neugebauer et al. 2022. **AH TODO: Add links/references to papers**
## Installation
We recommend installing everything in a virtual environment, though this is optional. To download this repository and create an appropriate environment, run:
```
git clone https://github.com/alvin-hsu/BELib_Processing
cd BELib_Processing
python3 -m venv ./venv
source ./venv/bin/activate
pip3 install -r requirements.txt
```

## Input Files
The two main inputs that need to be supplied (along with the actual `.fastq.gz` data files) are CSV files describing the library design and the experimental design. The library design file must contain the following columns (the corresponding file for the 12kChar library is provided on Git LFS as an example here):  **AH TODO: Add link to 12k library**
- `Index`: The ID of the library member. Must be the first column in the file.
- `gRNA`: The sequence of the spacer used (keep in mind an extra G must be used for U6 transcription; this is **not** included here).
- `original`: The original (i.e., unedited) sequence of the target site.

The experimental design file must contain the following columns (the corresponding file for Neugebauer et al. is provided on Git LFS as an example here):
- `name`: The name of the sample. Note that this should match what was used in the samplesheet for `bcl2fastq` demultiplexing!
- `sample_num`: The sample number. Note that this should match what was used in the samplesheet for `bcl2fastq` demultiplexing!
- `cell_type`: The cell type used in the experiment.
- `editor_name`: The name of the base editor being examined.
- `rep`: The ID of the experimental replicate.

## Configuring the settings
For simplicity, settings used for processing are kept in `_config.py`. The settings provided here are the same as those used in Neugebauer et al. The only settings that should need checking each time are:
- `N_WORKERS`: Set this to the number of CPUs you want to use for alignment. If set to 1 or lower, reads are aligned serially (*very* slow!).
- `TARGET_START`, `TARGET_END`, `TARGET_REVCOMP`: Set these to the indices in the reverse read for where the target site begins and ends, as well as to whether the target site is reverse complemented (relative to the read).
- `GRNA_START`, `GRNA_END`, `GRNA_REVCOMP`: Set these to the indices in the forward read for where to look for the gRNA spacer, as well as to whether the gRNA is reverse complemented (relative to the read).
- `PROTOSPACER_OFFSET`: Set this to the offset in the (aligned) target site to describe where the protospacer is located.
- `BASE_EDITS`: Set this to dinucleotides describing the base edits you are interested in analyzing. For example, `AG` refers to an A-to-G edit.
- `WINDOW_START`, `WINDOW_END`: Set this to the positions (relative to the protospacer, such that position 1 is the start of the protospacer) that should be examined for base editing activity.

## Data processing
First, genotype the reads by running `1_align.py`. This script assigns reads from **one** `.fastq.gz` file to library members (we recommend running this script in parallel for each `.fastq.gz` file on a cluster/AWS). The usage is as follows:  **AH TODO: Add usage**
