# Usage
The main Python script VariantFastaConda.py aims at transforming a vcf file from the main scripts PairedEndVariantConda.py and VariantAlignmentConda.py into a multi-fasta file for downtream phylogenomic analyses, including unvariable sites from the single contig reference genome, as well as single nucleotide polymorphisms (SNPs) and/or small Insertions/Deletions (InDels) from each samples.
- The present main script VariantFastaConda.py corresponds to an adaptation with Python3 of the main script VCFtoPseudoGenome.py that Arnaud Felten developped with Python2.
- The main script VariantFastaConda.py and module genomic.py (version 20201006, October 2020) were prepared and tested with Python and dependencies below.
- The module genomic.py has to be with the present main script VariantFastaConda.py to launch it properly.
- The Conda environment PairedEndVariantCalling has to be prepared as presented below.
- The user can setup his own dependencies in his own bin.
- The input indexed reference and .vcf file must be preferably prepared with the main scripts PairedEndVariantConda.py and VariantAlignmentConda.py, respectively.
- The user can use as input his own .vcf file and indexed reference.
# Dependencies
The main script VariantFastaConda.py and module genomic.py (version 20201006) were prepared and tested with Conda packages below (Name/Version/Build/Channel).
- python/3.8.5/h1103e12_9_cpython/conda-forge
- biopython/1.78/py38h1e0a361_0/conda-forge
# Building of the Conda Environment PairedEndVariantCalling
## 1/ From available targeted Conda packages
```
conda activate
conda create -n PairedEndVariantCalling
conda activate PairedEndVariantCalling
conda search python
conda install -c conda-forge python=3.8.5=h1103e12_9_cpython
conda search biopython
conda install -c conda-forge biopython=1.78=py38h1e0a361_0
```
## 2/ From available updated Conda packages
```
conda activate
conda create -n PairedEndVariantCalling
conda activate PairedEndVariantCalling
conda install -c conda-forge python
conda update -c conda-forge python
conda install -c conda-forge biopython
conda update -c conda-forge biopython
```
# Launching of the script VariantFastaConda.py
## 1/ prepare a single command in a Bash script (bash_VariantFastaConda.sh)
```
#!/bin/bash
#SBATCH -p Research
#SBATCH -o %x.%N.%j.out
#SBATCH -e %x.%N.%j.err
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test-20201012
source /global/conda/bin/activate;conda activate PairedEndVariantCalling; \
python /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantFastaConda.py \
	-i /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/7_alignment/filtered.snps.indels.vcf \
	-ref /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/1_reference/Enteritidis_P125109.fasta \
	-o /global/bio/projets/GAMeR/Nicolas-Radomski/PairedEndVariant/VariantCalling/7_alignment/results.fasta \
	--NoINDELs
```
## 2/ run the Bash script bash_VariantFastaConda.sh with sbatch
```
sbatch bash_VariantFastaConda.sh
```
# Illustration
![Workflow](https://github.com/Nicolas-Radomski/VariantFastaConda/blob/main/illustration.png)
# References
- First version (i.e. GATK-SNP calling): Lee R.S., N. Radomski, J.F. Proulx, I. Levade, B.J. Shapiro, F. McIntosh, H. Soualhine, D. Menzies and M.A. Behr. Population genomics of Mycobacterium tuberculosis in the Inuit. 2015, Proceedings of the National Academy of Sciences of the United States of America, PNAS, 112(44): 13609-13614, doi: 10.1073/pnas.1507071112
- Second version (i.e. VarCall): Felten A., M. Vila Nova, K. Durimel, L. Guillier, M.Y. Mistou and N. Radomski. First gene-ontology enrichment analysis based on bacterial coregenome variants: insights into adaptations of Salmonella serovars to mammalian- and avian-hosts. 2017, BMC Microbiology, 17(222): 1-20, doi.org/10.1186/s12866-017-1132-1
- Third verion (i.e. iVarCall2): Vila Nova M, Durimel K., La K., Felten A., Bessières P., Mistou M.Y., Mariadassou M. and N. Radomski. Genetic and metabolic signatures of Salmonella enterica subsp. enterica associated with animal sources at the pangenomic scale. 2019, BMC Genomics, 20(1): 814, doi: 10.1186/s12864-019-6188-x
# Acknowledgment
My old colleagues Arnaud Felten and Ludovic Mallet with whom I learned a lot about Python
# Author
Nicolas Radomski
