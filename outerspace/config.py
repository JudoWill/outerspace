"""configuration for defining motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from tomlkit import comment
from tomlkit import document
from tomlkit import table
from tomlkit.toml_file import TOMLFile
from argparse import ArgumentParser


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

    def __str__(self):
        return f'Cfg({self.filename})'



 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
