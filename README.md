# RE2C

- RE2C is a free, open-source meta-analysis software tool for genome-wide association study analysis, designed to perform both basic and advanced meta-analytic methods in an efficient manner. RE2C allows users to test heterogeneous effect size between individual statistics by using an optimization technique, spectral decomposition. Our optimization procedure is fast and robust and takes only a day for 1 000 000 SNPs.

## Getting Started

In order to download `RE2C`, you should clone this repository via the command
```
git clone git://github.com/cuelee/RE2C.git
cd RE2C
```

## Updating RE2C
You can update to the newest version of `RE2C` using git. First, navigate to your RE2C/ directory, then run
```
git pull
```

## User Guide
- Usage 
Command line arguments: bash RE2C.bash

|Option|Description|
|---|---|
|--i [FILE]|Input file (Required)|
|--o [FILE]|Output file (defulat = output/out.txt)|
|--c [FILE] or --cor [FILE]|Correlation matrix|
|--h or --help|Print help|

- Input File Format 

	+ Each rows represent a SNPs. 
		+ The first column of each rows is for rsID, 
and the following columns are pairs of effect size and its standard error of Nth study. (If we meta-analyze 5 summary statistics then each line must have 11 columns)

	+ Correlation matrix is N x N symmetric matrix. 
		+ An evaluated correlation matrix can be specified by using --cor option in RE2C. When --cor is unused, the code itself assumes an identity matrix of N rows and columns.

Example (Input File)
```
rsAAAAAA study1beta study1stderr study2beta study2stderr study3beta study3stderr
rsBBBBBB study1beta study1stderr study2beta study2stderr study3beta study3stderr
rsCCCCCC study1beta study1stderr study2beta study2stderr study3beta study3stderr
```
- Example (Correlation Matrix)

|C1|C2|C3|C4|
|---|---|---|---|
|1|a|b|c|
|a|1|d|e|
|b|d|1|h|
|c|e|h|1|

- Example Run
```
bash RE2C.sh --input ./example/example_input.txt --output ./out --cor ./example/example_cor.txt
```

- Output File

|Col. Num.|Col. name|Description|
|---|---|---|
|1|rsID|SNP rsID|
|2|Nstudy|Number of studies used in meta-analysis|
|3|LSp|Lin-Sullivan’s P-value (FE)|
|4|RE2Cs1|RE2C statistic mean effect part|
|5|RE2Cs2|RE2C statistic heterogeneity part|
|6|RE2Cp|RE2C p-value|
|7|RE2Cp*|RE2C* p-value|

## RE2C generates two statistical significances: RE2Cp and RE2Cp*. which one should I use for my data?

The original RE2C model (RE2Cp) was developed to complement the FE model, the golden standard in meta-analysis. For example, Let assume that we ran association tests on a SNP using 1. the conventional fixed effect meta-analysis, and 2. RE2C. If the p-value of the fixed effect method (FEp) is 0.03 and the p-value of the RE2C model (RE2Cp) is 0.01, the statistical significance of this SNP will be min(0.03, 0.01) * 2 = 0.02.

We later discovered that some might just want to use the statistical significance of RE2C alone. In this case, you can use RE2Cp* (not RE2Cp).

## Citation

If you use our software or refer to the RE2C method, please cite.
Lee CH, Eskin E, Han B. Increasing the power of meta-analysis of genome-wide association studies to detect heterogeneous effects. Bioinformatics. 2017;33(14):i379–88. 

## Support

Issues with RE2C? Email cuelee@snu.ac.kr

## License 

gnu license 3.0

## Authors

Cue Hyunkyu Lee ( Seoul National University )

## Funding Information

C.H.L B.H are supported by a grant of the Korean Health Technology R&D project, Ministry of Health & Welfare,
Republic of Korea (Hl14C1731).

