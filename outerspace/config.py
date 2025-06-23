"""configuration for defining motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from tomlkit import comment
from tomlkit import document
from tomlkit import table
from tomlkit.toml_file import TOMLFile
from argparse import ArgumentParser
from outerspace.pattern import Pattern


class Cfg:
    """configuration for defining motifs"""

    def __init__(self, filename=None):
        self.filename = filename

    @staticmethod
    def _generate_from_parser(parser: ArgumentParser, return_doc: bool = False):
        """Generate a TOML configuration file with default values from argparse parser.
        
        Args:
            parser: The ArgumentParser instance to generate config from
            return_doc: If True, return the TOML document object instead of string
            
        Returns:
            str | document: The TOML configuration file as a string or document object
        """
        doc = document()
        
        # Process each subcommand
        for cmd_name, subparser in parser._subparsers._group_actions[0].choices.items():
            section = table()
            
            # Extract defaults from subparser
            for action in subparser._actions:
                if action.dest != 'help':  # Skip help action
                    # Add help text and type as comment
                    help_text = action.help if action.help else "No description available"
                    type_info = f"Type: {action.type.__name__}" if action.type else "Type: str"
                    section.add(comment(f"\n# {help_text}\n# {type_info}"))
                    
                    if action.default is not None and action.default != '==SUPPRESS==':
                        section[action.dest] = action.default
                    elif action.default is None:
                        if not action.option_strings:  # Positional argument
                            section[action.dest] = 'positional'
                        else:  # Keyword argument
                            section[action.dest] = 'required' if action.required else 'optional'
            
            # Add section to document if it has any values
            if section:
                doc[cmd_name] = section
        
        return doc if return_doc else doc.as_string()

    def read_file(self):
        """Read the file specified"""
        return TOMLFile(self.filename).read() if self.filename is not None else None

    def write_file(self):
        """Write a default configuration file"""
        doc = self.get_doc_default()
        if self.filename is not None:
            TOMLFile(self.filename).write(doc)
        return doc

    @staticmethod
    def get_doc_default():
        """Getting default config document object"""

        from outerspace.cli.main import Cli

        return Cfg._generate_from_parser(Cli._init_parser(), return_doc=True)

    @staticmethod
    def parse_patterns_from_config(config_data, global_patterns=None):
        """Parse Pattern objects from config data.
        
        Args:
            config_data: Dictionary containing config data with pattern configuration
            global_patterns: Optional dictionary of global patterns by name
            
        Returns:
            list: List of Pattern objects
            
        Raises:
            ValueError: If pattern configuration is invalid
        """
        patterns = []
        
        # Check for pattern_names in config (new format)
        if 'pattern_names' in config_data and global_patterns:
            for pattern_name in config_data['pattern_names']:
                if pattern_name not in global_patterns:
                    raise ValueError(f"Pattern '{pattern_name}' not found in global patterns")
                pattern_config = global_patterns[pattern_name]
                patterns.append(Cfg._create_pattern_from_config(pattern_config))
            return patterns
        
        # Check for use_all_patterns flag
        if config_data.get('use_all_patterns', False) and global_patterns:
            for pattern_config in global_patterns.values():
                patterns.append(Cfg._create_pattern_from_config(pattern_config))
            return patterns
        
        # Check for inline patterns (old format)
        if 'patterns' in config_data:
            for i, pattern_config in enumerate(config_data['patterns']):
                patterns.append(Cfg._create_pattern_from_config(pattern_config, i))
            return patterns
        
        # No patterns found
        return patterns

    @staticmethod
    def _create_pattern_from_config(pattern_config, index=None):
        """Create a Pattern object from config data.
        
        Args:
            pattern_config: Dictionary containing pattern configuration
            index: Optional index for error reporting
            
        Returns:
            Pattern: Pattern object
            
        Raises:
            ValueError: If pattern configuration is invalid
        """
        try:
            # Validate required fields
            required_fields = ['reg_expr', 'read', 'orientation', 'multiple']
            for field in required_fields:
                if field not in pattern_config:
                    idx_str = f" {index}" if index is not None else ""
                    raise ValueError(f"Pattern{idx_str}: Missing required field '{field}'")
            
            # Create Pattern object
            pattern = Pattern(
                reg_expr=pattern_config['reg_expr'],
                read=pattern_config['read'],
                orientation=pattern_config['orientation'],
                multiple=pattern_config['multiple']
            )
            return pattern
            
        except Exception as e:
            idx_str = f" {index}" if index is not None else ""
            raise ValueError(f"Pattern{idx_str}: {str(e)}")

    @staticmethod
    def parse_global_patterns(toml_doc):
        """Parse global patterns from TOML document.
        
        Args:
            toml_doc: TOML document object
            
        Returns:
            dict: Dictionary of patterns by name
        """
        global_patterns = {}
        
        if 'patterns' in toml_doc:
            for pattern_config in toml_doc['patterns']:
                if 'name' not in pattern_config:
                    raise ValueError("Global patterns must have a 'name' field")
                name = pattern_config['name']
                global_patterns[name] = pattern_config
        
        return global_patterns

    def __str__(self):
        return f'Cfg({self.filename})'


 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
