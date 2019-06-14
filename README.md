# SELFISH
SELFISH (Discovery of Differential Chromatin Interactions via a Self-Similarity Measure) is a tool for finding differential chromatin interactions
between two Hi-C contact maps. It uses self-similarity to model interactions 
in a robust way. For more information read the full 
paper: <a href="https://www.biorxiv.org/content/10.1101/540708v1?rss=1" target="_blank">**Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure**</a>. 

## Installation
#### TODO

## Usage
`$ selfish -f1 NPC_chr4_100kb.txt -f2 ES_chr4_100kb.txt -o ./ -r 100kb -b normals.KRNorm` 

## Parameters
| Short | Long | Meaning |
|---|---|---|
|_Required Parameters_| | |
| **-f1** | **--file1** | Location of contact map 1. (See below for format.) |
| **-f2** | **--file2** | Location of contact map 2. (See below for format.) |
| **-r** | **--resolution** | Resolution of the provided contact maps. |
| **-o** | **--outdir** | Directory to save the output. |
| _Optional Parameters_ | | |
| **-b** | **--biases** | Location of biases file. (See below for format.) |
| **-sz** | **--sigmaZero** | Sigma0 parameter for Selfish. Default is experimentally chosen for 5Kb resolution.|
| **-i** | **--iterations** | Iteration count parameter for Selfish. Default is experimentally chosen for 5Kb resolution.|
| **-v** | **--verbose** | Whether the program prints its progress while running. Default is True. |
| **-p** | **--plot** | Whether the program plots its results. Default is False. |
| **-lm** | **--lowmem** | Use float32 instead of float64. Uses less memory at the cost of precision. Default is False. |
| **-V** | **--version** | Shows the version of the tool. |

### Input Formats
#### Contact Maps
Contact maps need to have the following format. They must not have a header. 
Values must be separated by either a space, a tab, or a comma.
Midpoints must either be all divisible to resolution, or be already divided. (ie. if the resolution is 5.000, midpoint can be 10 or 50.000)

| Midpoint 1 | Midpoint 2 | Contact Count |
|---|---|---|
| 5000 | 65000 | 438 |
| 5000 | 85000 | 12 |
| ... | ... | ... |

#### Bias File
Bias file need to have the following format.
Bias file must use the same midpoint format as the contact maps.
Bias file must not have a header.

| Midpoint | Bias |
|---|---|
| 5000 | 1.012 |
| 65000 | 1.158 |
| ... | ... |

### Output
Output of Selfish is a matrix of p-values indicating the probability of differential conformation (Smaller values mean more significant.). X and Y coordinates indicate the bin midpoints.