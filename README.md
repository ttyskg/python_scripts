# python_script
Personal small python scripts

## `assemblingseq`

Create DNA sequences for assembling, such as for Gibson assembly,
polymerase cycling assembly (PCA), and so on. The created sequences will be
saved as a Text file

This script requres two args.

1. Number of DNA sequences to create
1. Sequence length

__Example:__
```bash
python aassemblingseq.py 5 20
```

__Requirements:__

* biopython
* [UNAFold](http://www.unafold.org/) (only require `hybrid-ss-min`)


## `randomseq`

Create DNA sequences that meet the condition. The created sequences will be
saved as a CSV file.  
This script requires four args.  

1. Number of DNA sequences to create
1. Sequence length
1. Minimum GC content
1. Maximum GC content

__Example:__
```bash
randomseq.py 3 20 45 70
```
