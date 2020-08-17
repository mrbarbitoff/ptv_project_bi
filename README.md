# Harnessing population-specific protein-truncating variants to improve the annotation of loss-of-function alleles

This repository contains all data and code pertinent to the analysis presented in the paper. This project is dedicated to the identification of population-specific PTVs (psPTVs) and building a predictive model to enhance annotation of low-confidence loss-of-function (LoF) alleles. 

psPTVs are strongly enriched for low-confidence putative LoF (pLoF) variants. An example of a notable psPTV that is a misannotated low-conf pLoF is a rs139297920 in the *PAX3* gene:

![alt text](https://github.com/mrbarbitoff/ptv_project_bi/blob/master/pax3.png)

`psPTV_project.Rmd` markdown file represents the body of the analysis and contains the complete explanation of all data analysis steps. All files included in this repository are used in this Markdown file.

`LoFfeR` folder contains the source code and reference materials of the **LoFfeR** toolkit that was developed to enhance the annotation of low-confidence pLoF alleles. See below for usage details.

`shet_estimation` directory contains scripts used for s-het inference and calculation of the PTV count distribution likelihood and the observation likelihood score. `IG_params.py` fits the inverse Gaussian distribution to the data; `estimate_S.py` estimates the best-fit s-het value, `CPL.py` calculates the likelihood of a PTV count distribution and conducts sampling.

`cassa_table.csv` is the main supporting table of the Cassa *et al.,* 2017 paper.

`clv*tsv` are tab-separated data files containing the validation dataset of benign and confident pathogenic ClinVar variants. These data files are used during model training and evaluation.

`covered_genes.tsv` is the file listing coordinates of the genes that have at least 30x mean coverage of coding regions in gnomAD.

`full_data_constr_pext.tsv` is the main tab-separated data files containing the final annotated set of gnomAD PTVs used throughout the analysis. 

`full_corr_wshet.tsv`, `full_out_boot_stats.csv`, `comprehensive_gene_annotation.tsv` are the tab-seprated files with the aggregated gene-level PTV counts. `comprehensive_gene_annotation.tsv` is the final file containing all necessary annotations.

`merged_counts.tsv` is the GENCODE v19-derived tab-separated file containing per-gene isoform and exon counts.

`genes_inheritance.txt` is a simplified OMIM-derived data file listing the disease associations of the genes.

## LoFfeR usage

To run LoFfeR, a tool to predict low-confidence pLoF variants in the VCF file, you will need R v.3.6+ with the randomForest package and Python 3.6+ with pandas, numpy packages installed.

### Installation

To install LoFfeR, simply clone this repository and grant execution permissions to the main executable file:

```
git clone https://github.com/mrbarbitoff/ptv_project_bi
cd ptv_project_bi
chmod +x ./LoFfeR/LoFfeR
```

### Running LoFfeR

To run annotation of low-confidence pLoF variants, please use the following command:

```
./LoFfeR/LoFfeR <VCF>
```

By default, LoFfeR expects a VEP-annotated VCF file to be passed as the first and only argument. The file should be annotated using the LOFTEE plugin (LoF), as only variants passing LOFTEE filter will be considered. LOFTEE fields are expected to be the trailing ones in each annotation. 

**Currently**, LoFfeR only supports GRCh37 (b37) VCF files. We are working to add the support for the GRCh38 reference as soon as possible.

LoFfeR results are written to the `results.tsv` tab-separated file.
