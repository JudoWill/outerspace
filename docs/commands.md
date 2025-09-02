# Commands

Below are the commands available in `outerspace`:

### `outerspace findseq`
Extracts sequences from various file formats based on configuration patterns. Features:
- Support for FASTQ, FASTA, SAM, and BAM files
- Single file or paired-end processing
- Global pattern configuration system
- Auto-detection of file formats
- Progress tracking and error handling

```bash
usage: outerspace findseq [-h] config [-1 READ1_FILENAME] [-2 READ2_FILENAME] [-o OUTPUT_FILENAME]
                 [--region REGION] [--fetch {mapped,unmapped,all}]
                 [--long-format] [--matches-only]

Extract sequences from files based on configuration patterns

positional arguments:
  config                Configuration file with search patterns

options:
  -h, --help            show this help message and exit
  -1 READ1_FILENAME, --read1_filename READ1_FILENAME
                        Input file for read 1 (FASTQ, FASTA, SAM, BAM) or single read file
  -2 READ2_FILENAME, --read2_filename READ2_FILENAME
                        Input file for read 2 (FASTQ, FASTA, SAM, BAM) for paired reads
  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
                        Output CSV file name
  --region REGION       SAM/BAM region specification (e.g., "chr1:1-1000")
  --fetch {mapped,unmapped,all}
                        SAM/BAM fetch mode
  --long-format         Output in long format (one row per pattern match instead of one row per read)
  --matches-only        Only output reads that have at least one pattern match
```

### `outerspace collapse`
Corrects barcodes in CSV files using UMI-tools clustering. Features:
- Supports multiple barcode columns
- Configurable mismatch tolerance
- Multiple clustering methods (cluster, adjacency, directional)
- Row limiting for testing
- Detailed metrics output

```bash
usage: outerspace collapse [-h] (--input-file INPUT_FILE | --input-dir INPUT_DIR)
                 (--output-file OUTPUT_FILE | --output-dir OUTPUT_DIR)
                 --columns COLUMNS [--mismatches MISMATCHES] [--sep SEP]
                 [--row-limit ROW_LIMIT] [--method {cluster,adjacency,directional}]
                 [--metrics METRICS] [--config CONFIG]

Correct barcodes in CSV files using UMI-tools clustering

options:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        Input CSV file to process
  --input-dir INPUT_DIR
                        Input directory containing CSV files to process
  --output-file OUTPUT_FILE
                        Output CSV file for corrected barcodes
  --output-dir OUTPUT_DIR
                        Output directory for corrected CSV files
  --columns COLUMNS     Column(s) containing barcodes to correct. Can be a single column or comma-separated list
  --mismatches MISMATCHES
                        Number of mismatches allowed for clustering (default: 2)
  --sep SEP             CSV separator (default: ',')
  --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  --method {cluster,adjacency,directional}
                        Clustering method to use (default: directional)
  --metrics METRICS     Output YAML file for metrics
  --config CONFIG       TOML configuration file containing command settings
```

### `outerspace count`
Counts unique barcodes per key value in CSV files. Features:
- Barcode and key column specification
- `--allowed-list` filtering
- Downsampling capability
- Detailed output with barcode lists
- Gini coefficient calculation for both barcodes and keys
- Metrics output in YAML format

```bash
usage: outerspace count [-h] (--input-file INPUT_FILE | --input-dir INPUT_DIR)
               (--output-file OUTPUT_FILE | --output-dir OUTPUT_DIR)
               --barcode-column BARCODE_COLUMN --key-column KEY_COLUMN
               [--sep SEP] [--row-limit ROW_LIMIT] [--allowed-list ALLOWED_LIST]
               [--detailed] [--downsample DOWNSAMPLE] [--random-seed RANDOM_SEED]
               [--metrics METRICS] [--config CONFIG]

Count unique barcodes per key value in CSV files

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file INPUT_FILE
                        Input CSV file to process
  -d INPUT_DIR, --input-dir INPUT_DIR
                        Input directory containing CSV files to process
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output CSV file for barcode counts
  -D OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for barcode counts
  -b BARCODE_COLUMN, --barcode-column BARCODE_COLUMN
                        Column containing barcodes
  -k KEY_COLUMN, --key-column KEY_COLUMN
                        Column to group by
  -s SEP, --sep SEP     CSV separator (default: ',')
  -l ROW_LIMIT, --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  -a ALLOWED_LIST, --allowed-list ALLOWED_LIST
                        Text file containing allowed keys (one per line)
  --detailed            Include barcode lists in output
  --downsample DOWNSAMPLE
                        Randomly sample reads with probability between 0 and 1
  --random-seed RANDOM_SEED
                        Random seed for downsampling
  -m METRICS, --metrics METRICS
                        Output YAML file for metrics
  -c CONFIG, --config CONFIG
                        YAML configuration file for command
```

### `outerspace gini`
Calculates Gini coefficient from counts in a CSV column. Features:
- Support for pre-counted values
- Scale factor application for dealing with 'pre-normalized' counts
- `--allowed-list` filtering
- Detailed statistics output

```


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.