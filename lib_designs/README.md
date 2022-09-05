## Library design files
These libraries were designed by Mandana Arbab and Max Shen. They are formatted
to contain the columns `name`, `gRNA`, and `original`. Here, `original` refers
to the unedited target sequence in the library (mostly relevant for the `AtoG`
and `CtoT` libraries).

If you wish to use one of these libraries, the alignment parameters in
`_config.py` should be set according to the following table:
| Library Name        | `OFFSET` |`TARGET_LEN`|
|---------------------|:--------:|:----------:|
|`library_12kChar.csv`| 21       | 56         |
|`library_AtoG.csv`   | 10       | 56         |
|`library_CtoT.csv`   | 10       | 56         |
|`library_LibA.csv`   |  9       | 56         |
|`library_PAMvar.csv` | 13       | 56         |

If using the same primers as reported previously, the `TARGET_START` should be
`0`, `TARGET_REVCOMP` is `True`, `GRNA_START` should be `18`, `GRNA_END` should
be `42`, and `GRNA_REVCOMP` is `False` (for clarity, see library design in the
Supplementary Information of Shen, Arbab, et al. Nature (2018)).

