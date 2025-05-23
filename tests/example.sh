#STEPS
## 1.Test this shell script to ensure no errors: Run this file on command line utilizing: `sh example.sh`
## 2.Once verified it works create configs/grnaquery.cfg file unique to your regex expression 
## 3.Make a copy of `example.sh` file and assure the following:
### - rename new file
### - modify paths to your reads directory or read files and `results` directories
### - create directories in 'results' if do not exist (ie. results/extracted, results/collapsed, results/counted) 
## 4.Run created file on command line `sh newfile.sh`

#Term:
#Key= targeted pattern of interest/ sequence you are searching for


# Extract relevant regions from sequencing reads for each experiment 
# Relevant regions defined by .cfg file
findseq configs/grnaquery.cfg -1 data/409-4_S1_L002_R1_001.fastq.gz -2 data/409-4_S1_L002_R2_001.fastq.gz -o results/extracted/shuffle.csv
findseq configs/grnaquery.cfg -1 data/2-G1L9-M1_S9_L001_R1_001.fastq.gz -2 data/2-G1L9-M1_S9_L001_R2_001.fastq.gz -o results/extracted/M1-lib.csv
findseq configs/grnaquery.cfg -1 data/2-G1L9-M2_S12_L001_R1_001.fastq.gz -2 data/2-G1L9-M2_S12_L001_R2_001.fastq.gz -o results/extracted/M2-lib.csv



# Collapse all files in the extracted directory into the collapsed directory
# Collapsing barcodes, mismatch threshold is for barcode not key (ie. protospacer)
python ../original/collapse.py --mismatches 2 --columns "UMI_5prime,UMI_3prime" results/extracted/ results/collapsed/



# Count how many barcodes for each unique key (ie. protospacer)
# Ex. for each unique protospacer it will count how many instances of different barcodes it sees
# Metrics file created additionally: `counts.yaml` 
python ../original/count.py --barcode-column "UMI_5prime_UMI_3prime_corrected" --key-column protospacer results/collapsed/ results/counted/ --metrics results/counted/counts.yaml



# Additional Option for count: 
# Create allowed list of keys, provide this list and it will only count the item if in this list (Counting keys not barcodes)
# Utilize `--allowed` to provide list:
# `python ../original/count.py --allowed-list results/library_protospacers.txt --barcode-column "UMI_5prime_UMI_3prime_corrected" --key-column protospacer results/collapsed/ results/counted_allowed/ --metrics results/counted_allowed/counts.yaml`

# Note: This is done because we want exact match for the protospacer which can be one protospacer apart but do not want to include mismatch error in here and get back false positives bc sequence error occurs which can cause a 1nt change





## TODO: Add 'diff' script to find keys enriched in M1 vs M2
