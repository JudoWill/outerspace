# Outerspace CLI API Documentation

## Overview

The Outerspace CLI is implemented using a modular command-based architecture. Each command is implemented as a subclass of `BaseCommand`, providing a consistent interface for command-line operations while allowing for command-specific functionality.

## Main CLI Entry Point

The main CLI entry point is defined in `outerspace/cli/main.py`. The `Cli` class handles:

- Argument parsing setup
- Command registration
- Command instantiation and execution

The main workflow is:
1. Parse command-line arguments
2. Initialize the appropriate command class
3. Execute the command

## Base Command Class

The `BaseCommand` class (`outerspace/cli/commands/base.py`) serves as the foundation for all commands. It provides:

### Core Methods

- `_init_parser(subparsers)`: Abstract method that each command must implement to define its command-line arguments
- `run()`: Abstract method that each command must implement to define its execution logic
- `_load_config(config_file)`: Loads settings from a TOML configuration file
- `_merge_config_and_args(defaults)`: Merges command-line arguments, config file settings, and defaults
- `_chk_exists(filenames)`: Utility method to verify file existence

### Separation of Concerns

The BaseCommand class handles common functionality that should be consistent across all commands, while leaving command-specific logic to the subclasses:

#### BaseCommand Responsibilities:
1. **Configuration Management**
   - Loading TOML configuration files
   - Merging settings from different sources
   - Tracking explicitly set arguments
   - Providing default value handling

2. **File Validation**
   - Checking file existence
   - Validating file types (file vs directory)
   - Providing consistent error messages

#### Subclass Responsibilities:
1. **Command-Specific Arguments**
   - Defining command-specific argument parsers
   - Setting up help text and argument groups
   - Defining required vs optional arguments

2. **Command Logic**
   - Implementing the actual command functionality
   - Processing input files
   - Generating output
   - Handling command-specific errors

3. **Command-Specific Validation**
   - Validating command-specific requirements
   - Checking argument combinations
   - Implementing command-specific file handling

### Configuration Management

The base command implements a sophisticated configuration system that merges settings from multiple sources. The merging process is handled by `_merge_config_and_args()` and follows a strict priority order:

1. **Command-line Arguments** (Highest Priority)
   - Any argument explicitly set via command line takes precedence
   - The system tracks which arguments were explicitly set by comparing against defaults
   - Once an argument is marked as explicit, it cannot be overridden by config or defaults

2. **Configuration File Settings** (Medium Priority)
   - Settings from the TOML config file are applied next
   - Only applied to arguments that weren't explicitly set via command line
   - Command-specific section is automatically detected based on class name
   - Example: `FindSeqCommand` looks for `[findseq]` section

3. **Default Values** (Lowest Priority)
   - Default values are applied to any remaining unset arguments
   - These are typically defined in the command's `run()` method
   - Serves as a fallback when no other value is specified

#### Configuration Merging Example:
```python
# In command subclass
def run(self):
    # Define defaults
    defaults = {
        'input': None,
        'output': 'output.txt',
        'threads': 1
    }
    
    # Merge with config and args
    self._merge_config_and_args(defaults)
    
    # At this point, self.args will have values from:
    # 1. Command line if explicitly set
    # 2. Config file if not explicit
    # 3. Defaults if neither of the above
```

#### TOML Configuration Example:
```toml
[findseq]
output = "default_output.txt"
threads = 4

[collapse]
input = "default_input.txt"
```

In this example:
- If `--threads 8` is specified on command line, it will use 8
- If not specified on command line but `threads = 4` in config, it will use 4
- If neither is specified, it will use the default of 1

## Implementing New Commands

To implement a new command:

1. Create a new class that inherits from `BaseCommand`
2. Implement the required abstract methods:
   ```python
   def _init_parser(self, subparsers):
       parser = subparsers.add_parser('commandname',
           help='Description of command')
       # Add command-specific arguments
       return parser

   def run(self):
       # Implement command logic
   ```

3. Register the command in `main.py`:
   - Add import statement
   - Add to command list in `_init_parser()`
   - Add to `command_map` in `_init_command()`

### Example Command Implementation

```python
class ExampleCommand(BaseCommand):
    def _init_parser(self, subparsers):
        parser = subparsers.add_parser('example',
            help='Example command')
        parser.add_argument('--input',
            help='Input file')
        parser.add_argument('--output',
            help='Output file')
        return parser

    def run(self):
        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)
        
        # Merge config and args with defaults
        defaults = {
            'input': None,
            'output': None
        }
        self._merge_config_and_args(defaults)
        
        # Implement command logic
```

## Configuration Files

Commands can use TOML configuration files to provide default settings. The configuration system:

1. Looks for a section matching the command name (lowercase, without 'Command' suffix)
2. Merges these settings with command-line arguments
3. Respects command-line argument precedence

### Configuration Generation

The `Cfg` class in `outerspace/config.py` provides functionality to automatically generate TOML configuration files from the argument parser definitions. This ensures that configuration files stay in sync with the command-line interface. Key features:

1. **Automatic Generation**
   - Creates TOML sections for each command
   - Includes help text and type information as comments
   - Preserves default values from argument definitions
   - Handles both positional and optional arguments

2. **Usage Example**
   ```python
   from outerspace.config import Cfg
   
   # Generate default configuration
   cfg = Cfg("config.toml")
   cfg.write_file()  # Creates config.toml with defaults
   ```

3. **Generated TOML Format**
   ```toml
   [commandname]
   # Help text for argument1
   # Type: str
   argument1 = "default_value"
   
   # Help text for argument2
   # Type: int
   argument2 = 42
   ```

This automatic generation ensures that:
- Configuration files are always up-to-date with the CLI interface
- Help text and type information is preserved in the config file
- Default values are consistently applied
- The configuration format is standardized across all commands

Example TOML configuration:
```toml
[example]
input = "default_input.txt"
output = "default_output.txt"
```

## Error Handling

The base command provides built-in error handling for:
- Missing configuration files
- Non-existent files
- Invalid directories
- Missing required arguments

## Best Practices

When implementing new commands:

1. Always implement both required abstract methods
2. Use the configuration system for default values
3. Validate inputs using `_chk_exists()` for file operations
4. Provide clear help text for all arguments
5. Handle both single-file and batch processing cases when appropriate
6. Use the configuration merging system to respect user preferences
7. Keep command-specific logic in the subclass
8. Leverage base class functionality for common operations
9. Document any command-specific configuration options
10. Use the built-in error handling where possible 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
