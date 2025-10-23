# panaroo

[![Build Status](https://github.com/gtonkinhill/panaroo/workflows/panaroo-CI/badge.svg)](https://github.com/gtonkinhill/panaroo/actions)

An updated update to the pipeline for pangenome investigation. This branch ONLY works for finding EXACT protein variants, and should ONLY be run with parameters as follows: 

```{bash}

panaroo -i <your_gff_files> -o <out_dir>
--clean-mode strict
--threshold 1.0
--family_threshold 1.0
--len_dif_percent 1.0
--family_len_dif_percent 1.0
--refind-mode off
--merge_paralogs (optional, but how I use it)
--remove-invalid-genes
-t <threads>

```

## Documentation

Documentation for Panaroo can be found [here](https://gtonkinhill.github.io/panaroo)
