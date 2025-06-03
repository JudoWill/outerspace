# Commands

Below are a proposed set of commands for `outerspace` based on the usages described.

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
                 [--metrics METRICS]

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
               [--metrics METRICS]

Count unique barcodes per key value in CSV files

options:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        Input CSV file to process
  --input-dir INPUT_DIR
                        Input directory containing CSV files to process
  --output-file OUTPUT_FILE
                        Output CSV file for barcode counts
  --output-dir OUTPUT_DIR
                        Output directory for barcode counts
  --barcode-column BARCODE_COLUMN
                        Column containing barcodes
  --key-column KEY_COLUMN
                        Column to group by
  --sep SEP             CSV separator (default: ',')
  --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  --allowed-list ALLOWED_LIST
                        Text file containing allowed keys (one per line)
  --detailed            Include barcode lists in output (default: False)
  --downsample DOWNSAMPLE
                        Randomly sample reads with probability between 0 and 1
  --random-seed RANDOM_SEED
                        Random seed for downsampling
  --metrics METRICS     Output YAML file for metrics
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

Calculate Gini coefficient from counts in a CSV column

positional arguments:
  input_file            Input CSV file

options:
  -h, --help            show this help message and exit
  --column COLUMN       Column to calculate Gini coefficient from
  --count-column COUNT_COLUMN
                        Column containing pre-counted values
  --scale SCALE         Scale factor for normalized values (e.g., if normalized to mean=1)
  --sep SEP             CSV separator (default: ',')
  --allowed-list ALLOWED_LIST
                        Text file containing allowed values (one per line)
```

### `outerspace visualize`
Creates visualizations of barcode counts from CSV files. Features:
- Histogram generation
- Configurable bins and scales
- Customizable titles and labels
- Log scale support
- Multiple output formats
- Statistical summary overlay

```bash
usage: outerspace visualize [-h] input_dir output_dir [--sep SEP] [--bins BINS]
                   [--title-prefix TITLE_PREFIX] [--xlabel XLABEL]
                   [--ylabel YLABEL] [--log-scale] [--format FORMAT]

Visualize barcode counts from CSV files

positional arguments:
  input_dir             Input directory containing CSV files with barcode counts
  output_dir            Output directory for visualization plots

options:
  -h, --help            show this help message and exit
  --sep SEP             CSV separator (default: ',')
  --bins BINS           Number of histogram bins (default: 50)
  --title-prefix TITLE_PREFIX
                        Prefix for plot titles (default: filename)
  --xlabel XLABEL       X-axis label (default: Number of Unique Barcodes)
  --ylabel YLABEL       Y-axis label (default: Count)
  --log-scale           Use log scale for y-axis
  --format FORMAT       Output image format (default: png)
```

### `outerspace pipeline`
A wrapper to complete all of the steps in a single command. Features:
- Automated processing of multiple FASTQ files
- Configurable barcode correction parameters
- Optional metrics generation
- Progress tracking
- Error handling

```bash
usage: outerspace pipeline [-h] config_file --input-dir INPUT_DIR --output-dir OUTPUT_DIR
                 [--allowed-list ALLOWED_LIST] [--mismatches MISMATCHES]
                 [--method {cluster,adjacency,directional}]
                 [--barcode-columns BARCODE_COLUMNS] [--key-column KEY_COLUMN]
                 [--sep SEP] [--metrics]

Run the complete OUTERSPACE pipeline

positional arguments:
  config_file           Configuration file with search patterns

options:
  -h, --help            show this help message and exit
  --input-dir INPUT_DIR
                        Directory containing paired FASTQ read files
  --output-dir OUTPUT_DIR
                        Output directory for all pipeline results
  --allowed-list ALLOWED_LIST
                        Text file containing allowed keys (one per line)
  --mismatches MISMATCHES
                        Number of mismatches allowed for clustering (default: 2)
  --method {cluster,adjacency,directional}
                        Clustering method to use (default: directional)
  --barcode-columns BARCODE_COLUMNS
                        Column(s) containing barcodes to correct (default: UMI_5prime,UMI_3prime)
  --key-column KEY_COLUMN
                        Column to group by for counting (default: protospacer)
  --sep SEP             CSV separator (default: ',')
  --metrics             Generate metrics files for collapse and count steps
```
