# SELFISH
SELFISH (Discovery of Differential Chromatin Interactions via a Self-Similarity Measure) is a tool by Abbas Roayaei Ardakany, Ferhat Ay, and Stefano Lonardi. It is currently maintained by Tuvan Gezer (hgezer@lji.org).  
SELFISH is a tool for finding differential chromatin interactions
between two Hi-C contact maps. It uses self-similarity to model interactions 
in a robust way. For more information read the full 
paper: <a href="https://www.biorxiv.org/content/10.1101/540708v1?rss=1" target="_blank">**Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure**</a>. 
![DCI](/demo.png)
## Installation and usage
### PIP
```bash
pip3 install selfish-hic
selfish -f1 /.../map1.txt \
        -f2 /.../map2.txt \
	-ch 2 \
        -r 100kb -o ./output.npy
```
### Github
Make sure you have Python 3 installed, along with all the dependencies listed.
```bash
git clone https://github.com/ay-lab/selfish
./selfish/selfish/selfish.py -f1 /.../map1.txt \
                             -f2 /.../map2.txt \
	                     -ch 2 \
                             -r 100kb -o ./output.npy
```
### Nextflow
If you have any problem regarding dependencies or version mismatches, we recommend using <a href="https://www.nextflow.io/" target="_blank">Nextflow</a> with a container technology like <a href="https://www.docker.com/get-started" target="_blank">Docker</a> or <a href="https://singularity.lbl.gov/" target="_blank">Singularity</a>. These methods require Nextflow(Can be installed with a single command that doesn't require special permissions.), and the desired container technology to be available.
Program arguments are given to Nextflow with two dashes and the short format listed below.   
**Updating:** If Nextflow warns that your project is outdated, use `nextflow pull ay-lab/selfish` in order to update to latest version.
#### Install Nextflow
Nextflow works for Linux and OS X. Install it using **one** of the commands listed below. **Requires Java 8+**
```bash
wget -qO- https://get.nextflow.io | bash
OR
curl -s https://get.nextflow.io | bash
```

#### With Docker
```bash
./nextflow run ay-lab/selfish --f1="/.../map1.txt" \
                              --f2="/.../map2.txt" \
	                      --ch 2 \
                              --r=100kb -profile docker
```

#### With Singularity
```bash
./nextflow run ay-lab/selfish --f1="/.../map1.txt" \
                              --f2="/.../map2.txt" \
	                      --ch 2 \
                              --r=100kb -profile singularity
```

### Bioconda
Bioconda install isn't currently available.
### Dependencies
Selfish uses some python packages to accomplish its mission. These are the packages used by selfish:
1. numpy
2. pandas
3. matplotlib
4. seaborn
5. scipy
6. statsmodels
7. hic-straw
8. pathlib

## Parameters
| Short | Long | Meaning |
|---|---|---|
|_Required Parameters_| | |
| **-f1** | **--file1** | Location of contact map 1. (See below for format.) Not required for HiC-Pro input. |
| **-f2** | **--file2** | Location of contact map 2. (See below for format.) Not required for HiC-Pro input. |
| **-r** | **--resolution** | Resolution of the provided contact maps. |
| **-o** | **--outfile** | Name of the output file. |
| **-ch** | **--chromosome** | Specify which chromosome to run the program for. |
| **-m1** | **--matrix1** | Location of matrix file, **only** for HiC-Pro type input. |
| **-m2** | **--matrix2** | Location of matrix file, **only** for HiC-Pro type input. |
| **-bed1** | **--bed1** | Location of bed file, **only** for HiC-Pro type input. |
| **-bed2** | **--bed2** | Location of bed file, **only** for HiC-Pro type input. |
| _Optional Parameters_ | | |
| **-b1** | **--biases1** | Location of biases file for contact map 1. (See below for format.) |
| **-b2** | **--biases2** | Location of biases file for contact map 2. (See below for format.) |
| **-sz** | **--sigmaZero** | Sigma0 parameter for Selfish. Default is experimentally chosen for 5Kb resolution.|
| **-i** | **--iterations** | Iteration count parameter for Selfish. Default is experimentally chosen for 5Kb resolution.|
| **-v** | **--verbose** | Whether the program prints its progress while running. Default is True. |
| **-p** | **--plot** | Whether the program plots its results. Highly discouraged for high resolutions(<50kb) as it will take a lot of time to compute the plots. For high resolutions, we recommend using the output matrix and plotting small sections of it manually. Default is False. |
| **-lm** | **--lowmem** | Use float32 instead of float64. Uses less memory at the cost of precision. Default is False. |
| **-V** | **--version** | Shows the version of the tool. |

### Input Formats
SELFISH supports 3 different input formats. **Plain text**, **.hic**, **.bed/.matrix** pairs(HiC-Pro format).
#### Text Contact Maps
Contact maps need to have the following format. They must not have a header. 
Values must be separated by either a space, a tab, or a comma.

| Chromosome | Midpoint 1 | Chromosome | Midpoint 2 | Contact Count |
|---|---|---|---|---|
| chr1 | 5000 | chr1 | 65000 | 438 |
| chr1 | 5000 | chr1 | 85000 | 12 |
| ... | ... | ... | ... | ... |


#### Bed-Matrix pairs (HiC-Pro format)
**User must provide a chromosome with the -ch argument.**  
.bed and .matrix files must have the same name other than the extension.
 Either file name can be provided as an input for selfish and the program will search for the second file automatically.
 
 
#### HiC file
**User must provide a chromosome with the -ch argument.**  
Selfish uses juicer's *straw* tool to read .hic files.

#### Bias File
Bias file need to have the following format.
Bias file must use the same midpoint format as the contact maps.
Bias file must not have a header.

| Chromosome | Midpoint | Bias |
|---|---|---|
| chr1 | 5000 | 1.012 |
| 1 | 65000 | 1.158 |
| ... | ... | ... |

### Output
Output of Selfish is a matrix of p-values indicating the probability of differential conformation (Smaller values mean more significant.). 
X and Y coordinates indicate the bin midpoints.
File format of the output is a binary numpy file. It can be read by using Numpy as follows.  
```python
import numpy as np
matrix = np.load("/path/to/output/selfish.npy")
```
