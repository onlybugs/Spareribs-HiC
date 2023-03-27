# Spareribs-HiC: a One-stop, User-friendly, and Scalable Hi-C Data Analysis Pipeline


## Summary
Spareribs-HiC is implemented by snakemake, the most popular process writing language, which enables one-stop sequencing of downstream files `(.fq)` with various data preprocessing operations to obtain a classical representation of chromatin 3D structure data `(.cool/.mcool)`. In addition, Spareribs-HiC can automatically call embedded scripts `(.py/.R/.sh)` for chromatin compartment **(A/B)** analysis and topologically associating domain **(TAD)** analysis. More importantly, users can write their own personalized analysis scripts to be embedded in Spareribs-HiC to perform a variety of downstream analysis tasks in a synergistic manner.


## Dependency

Spareribs-HiC is written in snk with Snakemake framework.

The following versions are recommended when using iEnhance:

- Python
- R
- cooler
- bwa
- samtools
- OnTAD
- numpy
- matplotlib


## Data Preparation

### 1. fq data
For fastq data, we desire the input to be a high quality `.fq` file that have been quality controlled using `fastp` or `fastq+Trimmomaticz`. If your data are not **QC**, please use the above tools to QC and then into the folder `./in`.
## **_Note:_** All data should have the following naming format `{Non-duplicate_sample_name}_R*.fq`.


### 2. ref data
For *reference genome* data, we also desire it's ready to go. If you have nothing, go to the **UCSC** database and download the latest reference genome file `(.fa)` to `./ref` folder, then run the following command to build the index.
~~~shell
bwa index ref.fa # -p hg38
~~~


## Usage

### 1. Running
To execute the Spareribs-HiC, first use the following command to check the *DAG of the task execution*, and then execute the whole pipeline after checking that there are no errors.
~~~shell
snakemake -s sk.py --dag | dot -Tpdf > dag.pdf
~~~

Next, use the following command to complete the entire pipeline of Hi-C data analysis.
~~~shell
nohup snakemake -s sk.py -p -c 99 > hic.log &
~~~


### 2. Expanding
Spareribs-HiC is strongly scalable, and users can follow the steps below to personalize and expand the analysis.

1. Write the corresponding script code Collaborate with snakemake to complete I/O operations.

The following block of code can do the interaction between `.py` and `.snk`.
~~~python
inpfilen = snakemake.input[0]
...
np.savetxt(snakemake.output[0],cor_oe,delimiter="\t")
~~~

The following block of code can do the interaction between `.R` and `.snk`.
~~~R
TadRt = read.table(snakemake@input[[1]],sep = "\t",header = F)
...
outfilen = snakemake@output[[2]]
~~~

2. Put the debugged code into `./scripts` folder, then re-execute the .snk file.

~~~bash
snakemake -s sk.py -p -c 99
~~~

## **_Note:_** If you want to continue execution from a certain *RULE*, please refer to the official snakemake documentation https://snakemakecn.readthedocs.io/zh_CN/latest/.
