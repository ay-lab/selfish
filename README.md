# Selfish

Selfish (Discovery of Differential Chromatin Interactions via a Self-Similarity Measure) is a tool by Abbas Roayaei Ardakany, Ferhat Ay (ferhatay@lji.org), and Stefano Lonardi. This Python implementation is currently maintained by Tuvan Gezer (tgezer@sabanciuniv.edu). The original MATLAB implementation by Abbas (roayaei@gmail.com) can be found at https://github.com/ucrbioinfo/SELFISH.

SELFISH is a tool for finding differential chromatin interactions between two Hi-C contact maps. It uses self-similarity to model interactions in a robust way. For more information read the full
paper: <a href="https://academic.oup.com/bioinformatics/article/35/14/i145/5529135">Selfish: discovery of differential chromatin interactions via a self-similarity measure</a>.
![DCI](/demo.png)

## Installation and usage

### PIP

```bash
pip3 install selfish-hic
selfish -f1 /path/to/contact/map1.txt \
        -f2 /path/to/contact/map2.txt \
	-ch 2 \
        -r 100kb -o ./output.npy
```

### Github

Make sure you have Python 3 installed, along with all the dependencies listed.

```bash
git clone https://github.com/ay-lab/selfish
./selfish/selfish/selfish.py -f1 /path/to/contact/map1.txt \
                             -f2 /path/to/contact/map2.txt \
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
./nextflow run ay-lab/selfish --f1="/path/to/contact/map1.txt" \
                              --f2="/path/to/contact/map2.txt" \
	                      --ch=2 \
                              --r=100kb -profile docker
```

#### With Singularity

```bash
./nextflow run ay-lab/selfish --f1="/.../map1.txt" \
                              --f2="/.../map2.txt" \
	                      --ch=2 \
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

| Short                 | Long             | Meaning                                                                                                                                                                                                                                                            |
| --------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| _Required Parameters_ |                  |                                                                                                                                                                                                                                                                    |
| **-f1**               | **--file1**      | Location of contact map 1. (See below for format.) Not required for HiC-Pro input.                                                                                                                                                                                 |
| **-f2**               | **--file2**      | Location of contact map 2. (See below for format.) Not required for HiC-Pro input.                                                                                                                                                                                 |
| **-r**                | **--resolution** | Resolution of the provided contact maps.                                                                                                                                                                                                                           |
| **-o**                | **--outfile**    | Name of the output file.                                                                                                                                                                                                                                           |
| **-ch**               | **--chromosome** | Specify which chromosome to run the program for.                                                                                                                                                                                                                   |
| **-m1**               | **--matrix1**    | Location of matrix file, **only** for HiC-Pro type input.                                                                                                                                                                                                          |
| **-m2**               | **--matrix2**    | Location of matrix file, **only** for HiC-Pro type input.                                                                                                                                                                                                          |
| **-bed1**             | **--bed1**       | Location of bed file, **only** for HiC-Pro type input.                                                                                                                                                                                                             |
| **-bed2**             | **--bed2**       | Location of bed file, **only** for HiC-Pro type input.                                                                                                                                                                                                             |
| _Optional Parameters_ |                  |                                                                                                                                                                                                                                                                    |
| **-t**                | **--tsvout**     | If specified, outputs will be written as a TSV file. Specify the q-value threshold for which the results will be written to the file (i.e. -t 0.05)                                                                                                                |
| **-d**                | **--distanceFilter**    | Maximum distance between interacting loci.                                                                                                                                                                                          |
| **-c**                | **--changes**    | Name of the output file that has the log fold changes (log2((a+1)/(b+1))) between the inputs.                                                                                                                                                                                          |
| **-b1**               | **--biases1**    | Location of bias/normalization vector file for contact map 1. (See below for format.)                                                                                                                                                                                                 |
| **-b2**               | **--biases2**    | Location of bias/normalization vector file for contact map 2. (See below for format.)                                                                                                                                                                                                 |
| **-sz**               | **--sigmaZero**  | Sigma0 parameter for Selfish. Default is experimentally chosen for 5Kb resolution.                                                                                                                                                                                 |
| **-i**                | **--iterations** | Iteration count parameter for Selfish. Default is experimentally chosen for 5Kb resolution.                                                                                                                                                                        |
| **-v**                | **--verbose**    | Whether the program prints its progress while running. Default is True.                                                                                                                                                                                            |
| **-st**                | **--sparsityThreshold**     | The sparsity threshold selfish uses to filter differential interactions. Differential interactions falling in sparse regions are filtered out. The default is 0.7.                                                                                                                |
| **-ns**                | **--no-sparsityCheck**     | Do not use sparsity checking.                                                                                                                |
| **-nm**                | **--no-mutual**    | Include also the interactions that are zero in one of the contact maps. By defalut, selfish only considers interactions that are not zero in any of the contact maps.                                                                                                                                                                                            |
| **-norm**               | **--normalization**         | For .[m]cool or .hic files, you can specify what normalization selfish should use (-norm KR). For .[m]cool files, by default, selfish assumes that bias values are stored in the 'weight' column when this parameter is not specified. For .hic format, the default is 'KR' normalization if not specified.                                         |
| **-p**                | **--plot**       | Whether the program plots its results. Highly discouraged for high resolutions(<50kb) as it will take a lot of time to compute the plots. For high resolutions, we recommend using the output matrix and plotting small sections of it manually. Default is False. |
| **-lm**               | **--lowmem**     | Use float32 instead of float64. Uses less memory at the cost of precision. Default is False.                                                                                                                                                                       |
| **-V**                | **--version**    | Shows the version of the tool.                                                                                                                                                                                                                                     |

### Input Formats

SELFISH supports 3 different input formats. **Plain text**, **.hic**, **.cool**, **.bed/.matrix** pairs (HiC-Pro format).

#### Text Contact Maps

Contact maps need to have the following format. They must not have a header.
Values must be separated by either a space, a tab, or a comma.

| Chromosome | Midpoint 1 | Chromosome | Midpoint 2 | Contact Count |
| ---------- | ---------- | ---------- | ---------- | ------------- |
| chr1       | 5000       | chr1       | 65000      | 438           |
| chr1       | 5000       | chr1       | 85000      | 12            |
| ...        | ...        | ...        | ...        | ...           |

#### Bed-Matrix pairs (HiC-Pro format)

**User must provide a chromosome with the -ch argument.**
.bed and .matrix files must have the same name other than the extension.
Either file name can be provided as an input for selfish and the program will search for the second file automatically.

#### HiC file

**User must provide a chromosome with the -ch argument.**
Selfish uses juicer's _straw_ tool to read .hic files.

#### .cool file

**User must provide a chromosome with the -ch argument.**
Selfish uses _cooler_ to read .cool files.

#### Bias (normalization) File

Bias file need to have the following format.
Bias file must use the same midpoint format as the contact maps.
Bias file must not have a header.

| Chromosome | Midpoint | Bias  |
| ---------- | -------- | ----- |
| chr1       | 5000     | 1.012 |
| 1          | 65000    | 1.158 |
| ...        | ...      | ...   |

### Output

Output of Selfish is a matrix of q-values indicating the probability of differential conformation (Smaller values mean more significant.).
X and Y coordinates indicate the bin midpoints.  
Another optional output is the log fold changes file. It is simply produced by `log2((map1 + 1) / (map2 + 1))`  
File format of the outputs is a binary numpy file. It can be read by using Numpy as follows.

```python
import numpy as np
matrix = np.load("/path/to/output/selfish.npy")
```
## Citation
If you use Selfish in your work, please cite our <a href="https://academic.oup.com/bioinformatics/article/35/14/i145/5529135">Bionformatics paper</a>:

#### Abbas Roayaei Ardakany, Ferhat Ay, Stefano Lonardi, Selfish: discovery of differential chromatin interactions via a self-similarity measure, Bioinformatics, Volume 35, Issue 14, July 2019, Pages i145–i153, https://doi.org/10.1093/bioinformatics/btz362 
