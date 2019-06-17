# SELFISH
SELFISH (Discovery of Differential Chromatin Interactions via a Self-Similarity Measure) is a tool for finding differential chromatin interactions
between two Hi-C contact maps. It uses self-similarity to model interactions 
in a robust way. For more information read the full 
paper: <a href="https://www.biorxiv.org/content/10.1101/540708v1?rss=1" target="_blank">**Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure**</a>. 

## Installation and usage
### Github
Make sure you have python 3 installed, along with all the dependencies listed.
```bash
git clone https://github.com/tuvangezer/selfish
./selfish/selfish/selfish.py -f1 /.../map1.txt \
                             -f2 /.../map2.txt \
                             -r 100kb -o ./
```
### Nextflow
If you have any problem regarding dependancies or version mismatches, we recommend using <a href="https://www.nextflow.io/" target="_blank">Nextflow</a> with a container technology like <a href="https://www.docker.com/get-started" target="_blank">Docker</a> or <a href="https://singularity.lbl.gov/" target="_blank">Singularity</a>. These methods require Nextflow(Can be installed with a single command that doesn't require special permissions.), and the desired container technology to be available.
Program arguments are given to nextflow with two dashes and the short format listed below.   
**Updating:** If Nextflow warns that your project is outdated, use `nextflow pull tuvangezer/selfish` in order to update to latest version.
#### Install Nextflow
Nextflow works for Linux and OS X. Install it using **one** of the commands listed below. **Requires Java 8+**
```bash
wget -qO- https://get.nextflow.io | bash
OR
curl -s https://get.nextflow.io | bash
```

#### With Docker
```bash
./nextflow run tuvangezer/selfish --f1="/.../map1.txt" \
                                --f2="/.../map2.txt" \
                                --r=100kb -profile docker
```

#### With Singularity
```bash
./nextflow run tuvangezer/selfish --f1="/.../map1.txt" \
                                --f2="/.../map2.txt" \
                                --r=100kb -profile singularity
```

### Bioconda
Bioconda install isn't currently available.
### PIP
Pip install isn't currently available.
### Dependencies
Selfish uses some python packages to accomplish its mission. These are the packages used by selfish:
1. numpy
2. pandas
3. matplotlib
4. seaborn
5. scipy
6. statsmodels

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