# Differential_loops

#### Step1 generate chromatin loops and expected contacts

#### Step2 generate raw counts for chromatin loops and matrix for expected contacts, respectively

``` python GetCount_multiplesamples.py -I "mcool file path" -i "samples.tsv" -l "chromatin loop file" -O "output path" -o "output file" ```

usage: GetCount_multiplesamples.py [-h] -I INPUTPATH -i INPUT -l LOOPFILE -O OUTPATH -o OUTFILE

Build loop raw count and normalization factor matrices from mcool and expected files.

optional arguments:  
|:----:|:-----:|:----:|:------:|:------:|  
| -h | --help | show this help message and exit |  
| -I | INPUTPATH | --inputpath | INPUTPATH Directory containing the samples |  
| -i | INPUT | --input | INPUT Sample filename (tab-separated with columns: sample, mcool_path, resolution, expected_path) |  
| -l | LOOPFILE | --loopfile | LOOPFILE Loop file with BIN1_*/BIN2_* columns |  
|  -O | OUTPATH | --outpath | OUTPATH Output directory |  
|  -o | OUTFILE | --outfile | OUTFILE Output prefix (no extension) for the two result files |  

This code will generate EPloop_rawcounts.txt and EPloop_normalizationFactor.txt for Step3.
#### Step3 Run modified DESeq2 code to call differential chromatin loops

usage:
``` Rscript EPloop_DESeq2.R --Pathid "input and output file path" --Rawcounts "EPloop_rawcounts.txt" --Normalization "EPloop_normalizationFactor.txt" --Treatid LPS --Ctrlid Control --Treatnum 2 --Ctrlnum 2 --outfile LPS_vs_Control --log2FC 1 --padj 0.05 ```

optional arguments:  
|:----:|:------:|  
| --Pathid | Pathway of input and output files |  
| --Rawcounts | row-count matrix |  
| --Normalization | Size-factor file |  
| --Treatid | Treatment sampleid |  
| --Ctrlid | Control sampleid |  
| --Treatnum | Number of Treatment samples |  
| --Ctrlnum | Number of Control samples |  
| --outfile | Output filename |  
| --log2FC | log2FoldChange cutoff, default value is 1 |  
| --padj | padj cutoff, default value is 0.05 |  
