
# Extract relevant regions from sequencing reads for each experiment
findseq configs/grnaquery.cfg -1 data/409-4_S1_L002_R1_001.fastq.gz -2 data/409-4_S1_L002_R2_001.fastq.gz -o results/extracted/shuffle.csv
findseq configs/grnaquery.cfg -1 data/2-G1L9-M1_S9_L001_R1_001.fastq.gz -2 data/2-G1L9-M1_S9_L001_R2_001.fastq.gz -o results/extracted/M1-lib.csv
findseq configs/grnaquery.cfg -1 data/2-G1L9-M2_S12_L001_R1_001.fastq.gz -2 data/2-G1L9-M2_S12_L001_R2_001.fastq.gz -o results/extracted/M2-lib.csv

# Collapse all files in the extracted directory into the collapsed directory
python ../original/collapse.py --mismatches 2 --columns "UMI_5prime,UMI_3prime" results/extracted/ results/collapsed/

python ../original/count.py --barcode-column "UMI_5prime_UMI_3prime_corrected" --key-column protospacer results/collapsed/ results/counted/ --metrics results/counted/counts.yaml

python ../original/count.py --allowed-list results/library_protospacers.txt --barcode-column "UMI_5prime_UMI_3prime_corrected" --key-column protospacer results/collapsed/ results/counted_allowed/ --metrics results/counted_allowed/counts.yaml

## Diff script to find keys enriched in M1 vs M2
