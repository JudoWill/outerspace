# Commands

Below are the commands available in `outerspace`:

### `outerspace findseq`
Extracts sequences from FASTQ files based on configuration patterns. Features:
- Support for paired-end reads
- Configuration file for search patterns
- Batch processing of multiple read pairs
- Progress tracking with timing information
- Error handling for failed read pairs

```bash
usage: outerspace findseq [-h] config [-1 READ1_FILENAME] [-2 READ2_FILENAME] [-o OUTPUT_FILENAME]
                 [--fastqfiles FASTQFILES [FASTQFILES ...]] [--outdir OUTDIR]
                 [--read_regxlist READ_REGXLIST] [--read1_regxlist READ1_REGXLIST]
                 [--read2_regxlist READ2_REGXLIST]

Extract sequences from FASTQ files based on configuration patterns

positional arguments:
  config                Configuration file with search patterns

options:
  -h, --help            show this help message and exit
  -1 READ1_FILENAME, --read1_filename READ1_FILENAME
                        Zipped FASTQ file for read 1, or a single read
  -2 READ2_FILENAME, --read2_filename READ2_FILENAME
                        Zipped FASTQ file for read 2, or a single read
  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
                        Captured read file name output CSV
  --fastqfiles FASTQFILES [FASTQFILES ...]
                        Directory containing paired FASTQ read files
  --outdir OUTDIR       Output directory for processed files
  --read_regxlist READ_REGXLIST
                        Regular expression list for either read
  --read1_regxlist READ1_REGXLIST
                        Regular expression list for read 1
  --read2_regxlist READ2_REGXLIST
                        Regular expression list for read 2
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

```bash
usage: outerspace gini [-h] input_file --column COLUMN [--count-column COUNT_COLUMN]
              [--scale SCALE] [--sep SEP] [--allowed-list ALLOWED_LIST]
              [--config CONFIG]

Calculate Gini coefficient from counts in a CSV column

positional arguments:
  input_file            Input CSV file

options:
  -h, --help            show this help message and exit
  -c COLUMN, --column COLUMN
                        Column to calculate Gini coefficient from
  --count-column COUNT_COLUMN
                        Column containing pre-counted values
  --scale SCALE         Scale factor for normalized values (e.g., if normalized to mean=1)
  --sep SEP             CSV separator (default: ',')
  --allowed-list ALLOWED_LIST
                        Text file containing allowed values (one per line)
  --config CONFIG       TOML configuration file containing command settings
```

### `outerspace pipeline`
Runs the complete OUTERSPACE pipeline using Snakemake. Features:
- Automated processing of multiple FASTQ files
- Configurable barcode correction parameters
- Optional metrics generation
- Progress tracking
- Error handling

```bash
usage: outerspace pipeline [-h] config_file snakemake_config [--snakemake-args SNAKEMAKE_ARGS]

Run the complete OUTERSPACE pipeline using Snakemake

positional arguments:
  config_file           TOML configuration file with search patterns
  snakemake_config      YAML configuration file for Snakemake workflow

options:
  -h, --help            show this help message and exit
  --snakemake-args SNAKEMAKE_ARGS
                        Additional arguments to pass to Snakemake (e.g. --snakemake-args="--dry-run --cores 4")
```
