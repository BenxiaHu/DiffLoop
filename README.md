# Differential_loops

#### Step1 generate chromatin loops and expected contacts

#### Step2 generate raw counts for chromatin loops and matrix for expected contacts, respectively

python GetCount_multiplesamples.py -I "mcool file path" -i "samples.tsv" -l "chromatin loop file" -O "output path" -o "output file"

This code will generate EPloop_rawcounts.txt and EPloop_normalizationFactor.txt for Step3.
#### Step3 Run modified DESeq2 code to call differential chromatin loops

Rscript EPloop_DESeq2.R \
    --Pathid "input and output file path" \
    --Rawcounts "EPloop_rawcounts.txt" \
    --Normalization "EPloop_normalizationFactor.txt" \
    --Treatid LPS \
    --Ctrlid Control \
    --Treatnum 2 \
    --Ctrlnum 2 \
    --outfile LPS_vs_Control \
    --log2FC 1 \
    --padj 0.05