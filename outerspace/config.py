"""Configuration for defining motifs and managing TOML configuration files.

This module provides the Cfg class for handling configuration files in TOML format,
including parsing patterns, generating default configurations from argparse parsers,
and managing global pattern definitions.
"""
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SCB"

import logging
from typing import Any, Dict, List, Optional, Union
from tomlkit import comment, document, table
from tomlkit.toml_file import TOMLFile
from tomlkit.toml_document import TOMLDocument
from argparse import ArgumentParser

from outerspace.pattern import Pattern

# Set up logging
logger = logging.getLogger(__name__)


class Cfg:
    """Configuration manager for defining motifs and handling TOML configuration files.

    This class provides functionality to read, write, and parse TOML configuration
    files, generate default configurations from argparse parsers, and manage
    pattern definitions for motif analysis.
    """

    def __init__(self, filename: Optional[str] = None) -> None:
        """Initialize the configuration manager.

        Parameters
        ----------
        filename : Optional[str], default=None
            Path to the TOML configuration file to manage
        """
        self.filename = filename

    @staticmethod
    def _generate_from_parser(
        parser: ArgumentParser, return_doc: bool = False
    ) -> Union[str, TOMLDocument]:
        """Generate a TOML configuration file with default values from argparse parser.

        This method extracts argument information from an ArgumentParser instance
        and creates a corresponding TOML configuration structure with default values,
        help text, and type information.

        Parameters
        ----------
        parser : ArgumentParser
            The ArgumentParser instance to generate config from
        return_doc : bool, default=False
            If True, return the TOML document object instead of string

        Returns
        -------
        Union[str, TOMLDocument]
            The TOML configuration file as a string or document object

        Notes
        -----
        This method processes subcommands and their arguments, extracting defaults,
        help text, and type information to create a comprehensive configuration
        template.
        """
        doc = document()

        # Process each subcommand
        for cmd_name, subparser in parser._subparsers._group_actions[0].choices.items():
            section = table()

            # Extract defaults from subparser
            for action in subparser._actions:
                if action.dest != "help":  # Skip help action
                    # Add help text and type as comment
                    help_text = (
                        action.help if action.help else "No description available"
                    )
                    type_info = (
                        f"Type: {action.type.__name__}" if action.type else "Type: str"
                    )
                    section.add(comment(f"\n# {help_text}\n# {type_info}"))

                    if action.default is not None and action.default != "==SUPPRESS==":
                        section[action.dest] = action.default
                    elif action.default is None:
                        if not action.option_strings:  # Positional argument
                            section[action.dest] = "positional"
                        else:  # Keyword argument
                            section[action.dest] = (
                                "required" if action.required else "optional"
                            )

            # Add section to document if it has any values
            if section:
                doc[cmd_name] = section

        return doc if return_doc else doc.as_string()

    def read_file(self) -> Optional[TOMLDocument]:
        """Read the TOML configuration file.

        Returns
        -------
        Optional[TOMLDocument]
            The parsed TOML document, or None if no filename is set
        """
        if self.filename is None:
            logger.warning("No filename specified for reading configuration")
            return None
        return TOMLFile(self.filename).read()

    def write_file(self) -> TOMLDocument:
        """Write a default configuration file.

        Generates a default configuration document and writes it to the
        specified filename if one is set.

        Returns
        -------
        TOMLDocument
            The default configuration document that was written
        """
        doc = self.get_doc_default()
        if self.filename is not None:
            TOMLFile(self.filename).write(doc)
            logger.info(f"Default configuration written to {self.filename}")
        else:
            logger.warning("No filename specified for writing configuration")
        return doc

    @staticmethod
    def get_doc_default() -> TOMLDocument:
        """Get default configuration document object.

        Generates a default configuration by extracting argument information
        from the CLI parser.

        Returns
        -------
        TOMLDocument
            Default configuration document

        Notes
        -----
        This method imports the CLI module dynamically to avoid circular imports.
        """
        from outerspace.cli.main import Cli

        return Cfg._generate_from_parser(Cli._init_parser(), return_doc=True)

    @staticmethod
    def parse_patterns_from_config(
        config_data: Dict[str, Any], global_patterns: Optional[Dict[str, Any]] = None
    ) -> List[Pattern]:
        """Parse Pattern objects from configuration data.

        This method supports multiple configuration formats:
        - New format: Uses 'pattern_names' to reference global patterns
        - Legacy format: Uses 'use_all_patterns' flag or inline 'patterns'

        Parameters
        ----------
        config_data : Dict[str, Any]
            Dictionary containing configuration data with pattern configuration
        global_patterns : Optional[Dict[str, Any]], default=None
            Optional dictionary of global patterns by name

        Returns
        -------
        List[Pattern]
            List of Pattern objects parsed from the configuration

        Raises
        ------
        ValueError
            If pattern configuration is invalid or references non-existent patterns
        """
        patterns = []

        # Check for pattern_names in config (new format)
        if "pattern_names" in config_data and global_patterns:
            for pattern_name in config_data["pattern_names"]:
                if pattern_name not in global_patterns:
                    raise ValueError(
                        f"Pattern '{pattern_name}' not found in global patterns"
                    )
                pattern_config = global_patterns[pattern_name]
                patterns.append(Cfg._create_pattern_from_config(pattern_config))
            return patterns

        # Check for use_all_patterns flag
        if config_data.get("use_all_patterns", False) and global_patterns:
            for pattern_config in global_patterns.values():
                patterns.append(Cfg._create_pattern_from_config(pattern_config))
            return patterns

        # Check for inline patterns (old format)
        if "patterns" in config_data:
            for i, pattern_config in enumerate(config_data["patterns"]):
                patterns.append(Cfg._create_pattern_from_config(pattern_config, i))
            return patterns

        # No patterns found
        logger.debug("No patterns found in configuration data")
        return patterns

    @staticmethod
    def _create_pattern_from_config(
        pattern_config: Dict[str, Any], index: Optional[int] = None
    ) -> Pattern:
        """Create a Pattern object from configuration data.

        Parameters
        ----------
        pattern_config : Dict[str, Any]
            Dictionary containing pattern configuration
        index : Optional[int], default=None
            Optional index for error reporting

        Returns
        -------
        Pattern
            Pattern object created from configuration

        Raises
        ------
        ValueError
            If pattern configuration is invalid or missing required fields
        """
        try:
            # Validate required fields
            required_fields = ["reg_expr", "read", "orientation", "multiple"]
            for field in required_fields:
                if field not in pattern_config:
                    idx_str = f" {index}" if index is not None else ""
                    raise ValueError(
                        f"Pattern{idx_str}: Missing required field '{field}'"
                    )

            # Create Pattern object
            pattern = Pattern(
                reg_expr=pattern_config["reg_expr"],
                read=pattern_config["read"],
                orientation=pattern_config["orientation"],
                multiple=pattern_config["multiple"],
            )
            return pattern

        except Exception as e:
            idx_str = f" {index}" if index is not None else ""
            raise ValueError(f"Pattern{idx_str}: {str(e)}")

    @staticmethod
    def parse_global_patterns(toml_doc: TOMLDocument) -> Dict[str, Any]:
        """Parse global patterns from TOML document.

        Extracts pattern definitions from the 'patterns' section of a TOML
        document and returns them as a dictionary indexed by pattern name.

        Parameters
        ----------
        toml_doc : TOMLDocument
            TOML document object containing pattern definitions

        Returns
        -------
        Dict[str, Any]
            Dictionary of patterns by name

        Raises
        ------
        ValueError
            If global patterns are missing the required 'name' field
        """
        global_patterns = {}

        if "patterns" in toml_doc:
            for pattern_config in toml_doc["patterns"]:
                if "name" not in pattern_config:
                    raise ValueError("Global patterns must have a 'name' field")
                name = pattern_config["name"]
                global_patterns[name] = pattern_config

        return global_patterns

    def __str__(self) -> str:
        """Return string representation of the configuration manager.

        Returns
        -------
        str
            String representation showing the filename
        """
        return f"Cfg({self.filename})"


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
