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
usage: outerspace findseq [-h] [-1 READ1_FILENAME] [-2 READ2_FILENAME] [-o OUTPUT_FILENAME]
                 [--fastqfiles FASTQFILES [FASTQFILES ...]] [--outdir OUTDIR]
                 [config_filename]

Extract sequences from FASTQ files based on configuration patterns

positional arguments:
  config_filename       Configuration file with search patterns

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
usage: outerspace collapse [-h] --columns COLUMNS [--mismatches MISMATCHES] [--sep SEP]
                  [--row-limit ROW_LIMIT] [--method {cluster,adjacency,directional}]
                  input_dir output_dir

Correct barcodes in CSV files using UMI-tools clustering

positional arguments:
  input_dir             Input directory containing CSV files
  output_dir            Output directory for corrected CSV files

options:
  -h, --help            show this help message and exit
  --columns COLUMNS     Column(s) containing barcodes to correct. Can be a single column or comma-separated list
  --mismatches MISMATCHES
                        Number of mismatches allowed for clustering (default: 2)
  --sep SEP             CSV separator (default: ',')
  --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  --method {cluster,adjacency,directional}
                        Clustering method to use (default: directional)
```

### `count`
Counts unique barcodes per key value in CSV files. Features:
- Barcode and key column specification
- `--allowed-list` filtering
- Downsampling capability
- Detailed output with barcode lists
- Gini coefficient calculation for both barcodes and keys
- Metrics output in YAML format

```bash
usage: outerspace count [-h] --barcode-column BARCODE_COLUMN --key-column KEY_COLUMN
               [--sep SEP] [--row-limit ROW_LIMIT] [--allowed-list ALLOWED_LIST]
               [--detailed] [--downsample DOWNSAMPLE] [--random-seed RANDOM_SEED]
               [--metrics METRICS]
               input_dir output_dir

Count unique barcodes per key value in CSV files

positional arguments:
  input_dir             Input directory containing CSV files
  output_dir            Output directory for barcode counts

options:
  -h, --help            show this help message and exit
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

### `gini`
Calculates Gini coefficient from counts in a CSV column. Features:
- Support for pre-counted values
- Scale factor application for dealing with 'pre-normalized' counts
- `--allowed-list` filtering
- Detailed statistics output

```bash
usage: outerspace gini [-h] --column COLUMN [--count-column COUNT_COLUMN] [--scale SCALE]
              [--sep SEP] [--allowed-list ALLOWED_LIST]
              input_file

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
usage: outerspace visualize [-h] [--sep SEP] [--bins BINS] [--title-prefix TITLE_PREFIX]
                   [--xlabel XLABEL] [--ylabel YLABEL] [--log-scale]
                   [--format FORMAT]
                   input_dir output_dir

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

### `outspace pipeline`

A wrapper to complete all of the steps in a single command.
TBD