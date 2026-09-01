# The `DiffLoop` can identify Differential chromatin loops between 2 conditions.

#### Step1 generate chromatin loops and expected contacts

You can use loop callers and Cooltools to generate these files.


#### Step2 generate raw counts for chromatin loops and matrix for expected contacts, respectively

``` python GetCount_multiplesamples.py -I "mcool file path" -i "samples.tsv" -l "chromatin loop file" -O "output path" -o "output file" ```

usage:  
GetCount_multiplesamples.py [-h] -I INPUTPATH -i INPUT -l LOOPFILE -O OUTPATH -o OUTFILE

Build loop raw count and normalization factor matrices from mcool and expected files.

### Optional arguments

| Option              | Argument    | Description                                                                                            |
| :------------------ | :---------- | :----------------------------------------------------------------------------------------------------- |
| `-h`, `--help`      |             | Show the help message and exit                                                                         |
| `-I`, `--inputpath` | `INPUTPATH` | Directory containing the input files                                                                   |
| `-i`, `--input`     | `INPUT`     | Sample information file (tab-separated columns: `sample`, `mcool_path`, `resolution`, `expected_path`) |
| `-l`, `--loopfile`  | `LOOPFILE`  | Loop file containing `BIN1_*` and `BIN2_*` columns                                                     |
| `-O`, `--outpath`   | `OUTPATH`   | Output directory                                                                                       |
| `-o`, `--outfile`   | `OUTFILE`   | Output prefix (without extension) for the two result files                                             |


This code will generate EPloop_rawcounts.txt and EPloop_normalizationFactor.txt for Step3.


#### Step3 Run modified DESeq2 code to call differential chromatin loops

usage:  
``` Rscript EPloop_DESeq2.R --Pathid "input and output file path" --Rawcounts "EPloop_rawcounts.txt" --Normalization "EPloop_normalizationFactor.txt" --Treatid LPS --Ctrlid Control --Treatnum 2 --Ctrlnum 2 --outfile LPS_vs_Control --log2FC 1 --padj 0.05 ```

### Optional arguments

| Argument | Description | Type | Default |
|:---------|:------------|:----:|:-------:|
| `--Pathid` | Path to the input and output files | character | Required |
| `--Rawcounts` | Raw-count matrix file | character | Required |
| `--Normalization` | Normalization-factor file | character | Required |
| `--Treatid` | Treatment sample ID/prefix | character | Required |
| `--Ctrlid` | Control sample ID/prefix | character | Required |
| `--Treatnum` | Number of treatment replicates | integer | Required |
| `--Ctrlnum` | Number of control replicates | integer | Required |
| `--outfile` | Output file prefix | character | Required |
| `--log2FC` | Absolute log2 fold-change cutoff | numeric | `1` |
| `--padj` | Adjusted P-value cutoff | numeric | `0.05` |
