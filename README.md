# Base Editor Library Processing
Code for processing data from base editing screens performed with the libraries described in [Arbab,
Shen, et al., Cell (2020)](https://doi.org/10.1016/j.cell.2020.05.037). The code in this repository
was used to process data for Neugebauer et al., Nat. Biotechnol. (2022) **(in press)**.
## Installation
We recommend installing everything in a virtual environment, though this is optional. To download
this repository and create an appropriate environment, run the following code:
```
git clone https://github.com/alvin-hsu/BELib_Processing
cd BELib_Processing
python3 -m venv ./venv
source ./venv/bin/activate
pip3 install -r requirements.txt
```
If you want to use the library files provided in `lib_designs`, you will need to install
[git-lfs](https://git-lfs.github.com/) and run:
```
git lfs pull
```
## Input Files
The two main inputs that need to be supplied (along with the actual `.fastq.gz` data files) are CSV
files describing the library design and the experimental design. The library design file must
contain the following columns (library design files for libraries used in the Liu lab are available
in `lib_designs`):
- `Index`: The ID of the library member. Must be the first column in the file.
- `gRNA`: The sequence of the spacer used (keep in mind that an extra G might be
required for U6 transcription; this is **not** included here).
- `original`: The original (i.e., unedited) sequence of the target site.

The experimental design file must contain the following columns:
- `name`: The name of the sample. Note that this should match what was used in
the samplesheet for `bcl2fastq` demultiplexing!
- `sample_num`: The sample number. Note that this should match what was used in
the samplesheet for `bcl2fastq` demultiplexing!
- `cell_type`: The cell type used in the experiment.
- `editor_name`: The name of the base editor being examined.
- `rep`: The ID of the experimental replicate.

As an example, the corresponding file for Neugebauer et al. is provided on Git LFS as an example
[here](https://github.com/alvin-hsu/BELib_Processing/blob/main/exp_designs/exp_design_Neugebauer2022.csv)
(click "View raw").

If demultiplexing was done correctly, the required `.fastq.gz` files should have names like:
```
{name}_S{sample_num}_001.fastq.gz
```

## Configuring the settings
For simplicity, settings used for processing are kept in `_config.py`. The settings provided here
are the same as those used in Neugebauer et al. The only settings that should need checking each 
time are:
- `N_WORKERS`: Set this to the number of CPUs you want to use for alignment.
- `TARGET_START`, `TARGET_LEN`, `TARGET_REVCOMP`: Set these to the indices in the target read for
where the target site begins, the length of the target site, and whether the target site is reverse
complemented (relative to the read in which it occurs).
- `GRNA_START`, `GRNA_END`, `GRNA_REVCOMP`: Set these to the indices in the gRNA read for where to
look for the gRNA spacer, as well as to whether the gRNA is reverse complemented (relative to the 
read in which it occurs).
- `PS_OFFSET`: Set this to the offset in the (aligned) target site to describe where the protospacer
is located in the amplicon.
- `BASE_EDITS`: Set this to dinucleotides describing the base edits you are interested in analyzing.
For example, `AG` refers to an A•T-to-G•C edit.
- `WINDOW_START`, `WINDOW_END`: Set this to the positions (relative to the protospacer, such that
position 1 is the start of the protospacer) that should be examined for base editing activity.
- `ANALYSIS_NAMES`: Set this to a list of library member names that should be used for analyzing
activity. We need to set this because the library also contains members that were designed to look
at off-target editing.

## Data processing
First, genotype the reads by running `1_align.py`. This script assigns reads from **one** 
`.fastq.gz` file to library members (we recommend running this script in parallel for each 
`.fastq.gz` file on a cluster/AWS). The usage is as follows:
```
python 1_align.py [path_to_target_read] [path_to_grna_read] [path_to_lib_file]
```
Then, generate position-wise base editing summary files by running `2_process.py`. This script will
generate summaries for **all** library files. The output files `*_raw_stats.csv` are statistics 
generated from solely each genotype folder, while `*_corr_stats.csv` are statistics generated after
correcting for batch effects and Illumina sequencing error. The usage is as follows:
```
python 2_process.py [path_to_genotype_folders] [path_to_lib_design] [path_to_exp_design]
```

## Data Analysis
Finally, generate plots by running `3_analyze.py`, along with any other scripts desired. For
example, we have included the code to generate the comparison between TadCBEs with and without V106W
in `3_misc_Neugebauer.py`. The usage of `3_analyze.py` is identical to `2_process.py`:
```
python 3_analyze.py [path_to_genotype_folders] [path_to_lib_design] [path_to_exp_design]
```
