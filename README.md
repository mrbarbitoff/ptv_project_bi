# Harnessing population-specific protein-truncating variants to improve the annotation of loss-of-function alleles

This repository contains all data and code pertinent to the analysis presented in the paper. This project is dedicated to the identification of population-specific PTVs (psPTVs) and building a predictive model to enhance annotation of low-confidence loss-of-function (LoF) alleles. 

psPTVs are strongly enriched for low-confidence putative LoF (pLoF) variants. An example of a notable psPTV that is a misannotated low-conf pLoF is a rs139297920 in the *PAX3* gene:

![alt text](https://github.com/mrbarbitoff/ptv_project_bi/blob/master/pax3.png)

`psPTV_project.Rmd` markdown file represents the body of the analysis and contains the complete explanation of all data analysis steps and conclusions.

`LoFfeR` folder contains the source code and reference materials of the **LoFfeR** toolkit that was developed to enhance the annotation of low-confidence pLoF alleles. See below for usage details.

The other files and scripts included in this repository are:

## LoFfeR usage

To run LoFfeR, a tool to predict low-confidence pLoF variants in the VCF file, you will need R v.3.6+ (with the randomForest package), Python 3.6+ with pandas, numpy packages installed. please use the following command to run annotation:

```
git clone https://github.com/mrbarbitoff/ptv_project_bi
cd ptv_project_bi
chmod +x ./LoFfeR/LoFfeR
./LoFfeR/LoFfeR <VCF>
```

By default, LoFfeR expects a VEP-annotated VCF file to be passed as the first and only argument. The file should be annotated using the LOFTEE plugin (LoF), as only variants passing LOFTEE filter will be considered. LOFTEE fields are expected to be the trailing ones in each annotation. Currently, LoFfeR only supports GRCh37 (b37) VCF files. We are working to add the support for the GRCh38 reference as soon as possible.

LoFfeR results are written to the `results.tsv` tab-separated file.
