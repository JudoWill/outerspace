# Unified CLI Plan for GRNA Analysis Tools

## Overview

This document outlines the plan for combining the following commands into a single unified CLI interface:
- `findseq` (existing command)
- `collapse` (barcode correction)
- `count` (barcode counting)
- `gini` (Gini coefficient calculation)
- `visualize` (visualization)

The interface will follow a similar pattern to tools like samtools, where each command is a subcommand of the main tool.

## Command Structure

```bash
outerspace <command> [options]
```

Where `<command>` is one of:
- `find` - Extract sequences from fastq files
- `collapse` - Correct barcodes using UMI clustering
- `count` - Count unique barcodes per key
- `gini` - Calculate Gini coefficient
- `visualize` - Create visualizations

## Implementation Plan

### 1. Main CLI Class Structure

```python
class OuterspaceCLI:
    def __init__(self, args=None):
        self.parser = self._init_parser()
        self.args = self.parser.parse_args(args)
        self._prt_args()
        
    def run(self):
        """Route to appropriate command handler"""
        if not hasattr(self.args, 'command'):
            self.parser.print_help()
            sys_exit(1)
            
        command_handlers = {
            'find': self._run_find,
            'collapse': self._run_collapse,
            'count': self._run_count,
            'gini': self._run_gini,
            'visualize': self._run_visualize
        }
        
        handler = command_handlers.get(self.args.command)
        if handler:
            handler()
        else:
            self.parser.print_help()
            sys_exit(1)
```

### 2. Parser Structure

```python
def _init_parser(self):
    parser = ArgumentParser(
        prog='outerspace',
        description='GRNA analysis tools for sequence extraction and barcode analysis',
        epilog='Created by ThreeBlindMice - See how they run code'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Find command
    find_parser = subparsers.add_parser('find', help='Extract sequences from fastq files')
    find_parser.add_argument('config_filename', help='Configuration file with search patterns')
    find_parser.add_argument('-1', '--read1', help='Fastq file for read 1')
    find_parser.add_argument('-2', '--read2', help='Fastq file for read 2')
    find_parser.add_argument('-o', '--output', help='Output CSV file')
    find_parser.add_argument('--fastqfiles', nargs='*', help='Directory containing paired fastq files')
    find_parser.add_argument('--outdir', help='Output directory for processed files')
    
    # Collapse command
    collapse_parser = subparsers.add_parser('collapse', help='Correct barcodes using UMI clustering')
    collapse_parser.add_argument('input_dir', help='Input directory containing CSV files')
    collapse_parser.add_argument('output_dir', help='Output directory for corrected files')
    collapse_parser.add_argument('--columns', required=True, help='Column(s) containing barcodes')
    collapse_parser.add_argument('--mismatches', type=int, default=2, help='Number of mismatches allowed')
    collapse_parser.add_argument('--sep', default=',', help='CSV separator')
    collapse_parser.add_argument('--method', choices=['cluster', 'adjacency', 'directional'], 
                               default='directional', help='Clustering method')
    
    # Count command
    count_parser = subparsers.add_parser('count', help='Count unique barcodes per key')
    count_parser.add_argument('input_dir', help='Input directory containing CSV files')
    count_parser.add_argument('output_dir', help='Output directory for counts')
    count_parser.add_argument('--barcode-column', required=True, help='Column containing barcodes')
    count_parser.add_argument('--key-column', required=True, help='Column to group by')
    count_parser.add_argument('--allowed-list', help='File containing allowed keys')
    count_parser.add_argument('--downsample', type=float, help='Downsample probability')
    count_parser.add_argument('--metrics', help='Output YAML file for metrics')
    
    # Gini command
    gini_parser = subparsers.add_parser('gini', help='Calculate Gini coefficient')
    gini_parser.add_argument('input_file', help='Input CSV file')
    gini_parser.add_argument('--column', required=True, help='Column to calculate from')
    gini_parser.add_argument('--count-column', help='Column containing pre-counted values')
    gini_parser.add_argument('--scale', type=float, help='Scale factor for normalized values')
    gini_parser.add_argument('--allowed-list', help='File containing allowed values')
    
    # Visualize command
    visualize_parser = subparsers.add_parser('visualize', help='Create visualizations')
    visualize_parser.add_argument('input_dir', help='Input directory containing CSV files')
    visualize_parser.add_argument('output_dir', help='Output directory for plots')
    visualize_parser.add_argument('--bins', type=int, default=50, help='Number of histogram bins')
    visualize_parser.add_argument('--log-scale', action='store_true', help='Use log scale')
    visualize_parser.add_argument('--format', default='png', help='Output image format')
    
    return parser
```

### 3. Command Handlers

Each command will have its own handler method that implements the existing functionality:

```python
def _run_find(self):
    """Handle find command"""
    if self.args.config_filename is None:
        print('Please provide a config filename')
        sys_exit(1)
    # ... existing find command logic ...

def _run_collapse(self):
    """Handle collapse command"""
    # ... existing collapse command logic ...

def _run_count(self):
    """Handle count command"""
    # ... existing count command logic ...

def _run_gini(self):
    """Handle gini command"""
    # ... existing gini command logic ...

def _run_visualize(self):
    """Handle visualize command"""
    # ... existing visualize command logic ...
```

### 4. Entry Point

The main entry point will be updated to use the new CLI class:

```python
def main():
    """Main entry point for outerspace"""
    cli = OuterspaceCLI()
    cli.run()
```

## Usage Examples

```bash
# Extract sequences
outerspace find config.toml -1 read1.fastq -2 read2.fastq -o output.csv

# Correct barcodes
outerspace collapse input_dir output_dir --columns barcode --mismatches 2

# Count barcodes
outerspace count input_dir output_dir --barcode-column barcode --key-column sample

# Calculate Gini coefficient
outerspace gini input.csv --column counts

# Create visualizations
outerspace visualize input_dir output_dir --bins 50 --log-scale
```

## Benefits

1. **Unified Interface**: All tools are accessible through a single command
2. **Consistent Experience**: Common patterns and options across commands
3. **Better Discoverability**: Users can see all available commands with `outerspace --help`
4. **Simplified Installation**: Single package installation
5. **Shared Code**: Common functionality can be shared between commands

## Implementation Notes

1. Each command's functionality will remain largely unchanged
2. The existing code will be refactored into the new structure
3. Common utilities and helper functions can be shared
4. Error handling and logging will be consistent across commands
5. Documentation will be updated to reflect the new interface 