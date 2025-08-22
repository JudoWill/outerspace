"""Base class for all CLI commands.

This module provides the BaseCommand class that serves as the foundation for all
OUTERSPACE CLI commands. It handles common functionality including configuration
loading, argument parsing, pattern management, and file validation.
"""

import logging
import os
from argparse import ArgumentParser, Namespace
from pathlib import Path
from sys import exit as sys_exit
from typing import Any, Dict, Iterable, List, Optional, Set, Union
from itertools import islice

from tqdm import tqdm
from tomlkit import parse as toml_parse

from outerspace.config import Cfg
from outerspace.pattern import Pattern

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class BaseCommand:
    """Base class for all OUTERSPACE CLI commands.

    This class provides common functionality for all commands including:
    - Configuration file loading and parsing
    - Pattern management from TOML files
    - Argument merging (command line > config > defaults)
    - File existence validation
    - Common utility methods

    Subclasses must implement:
    - _init_parser(): Set up command-specific arguments
    - run(): Execute the command logic
    """

    def __init__(self, args: Optional[Namespace] = None) -> None:
        """Initialize the base command.

        Parameters
        ----------
        args : Optional[Namespace], default=None
            Parsed command-line arguments
        """
        self.args = args
        self._config: Optional[Dict[str, Any]] = None
        self._explicit_args: Set[str] = set()
        self._toml_doc: Optional[Any] = None
        self._global_patterns: Optional[Dict[str, Any]] = None

        logger.debug(f"Initialized {self.__class__.__name__}")

    def _init_parser(self, subparsers: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparsers : ArgumentParser
            Subparser group to add command arguments to

        Raises
        ------
        NotImplementedError
            This method must be implemented by subclasses
        """
        raise NotImplementedError("Each command must implement _init_parser")

    def _add_common_args(self, subparser: ArgumentParser) -> None:
        """Add common arguments to the command-specific parser.
        
        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """

        group = subparser.add_argument_group("Common Arguments")
        group.add_argument(
            "--config", '-c',
            type=str,
            help="Configuration file",
            default=None,
        )

        # Progress bar
        group.add_argument(
            "--progress-bar", '-p',
            action="store_true",
            help="Enable progress bar",
            default=False,
        )

        # Log file
        group.add_argument(
            "--log-file",
            type=str,
            help="Log file",
            default=None,
        )

        # Log level
        group.add_argument(
            "--log-level",
            type=str,
            help="Log level",
            default='WARNING',
        )

    def run(self) -> None:
        """Execute the command.

        Raises
        ------
        NotImplementedError
            This method must be implemented by subclasses
        """
        raise NotImplementedError("Each command must implement run")

    def _load_config(self, config_file: Union[str, Path]) -> Dict[str, Any]:
        """Load settings from TOML configuration file.

        This method loads a TOML configuration file, parses global patterns,
        and extracts command-specific settings based on the command class name.

        Parameters
        ----------
        config_file : Union[str, Path]
            Path to the TOML configuration file

        Returns
        -------
        Dict[str, Any]
            Command-specific configuration dictionary

        Raises
        ------
        ValueError
            If configuration file doesn't exist or is invalid
        """
        config_path = Path(config_file)

        if not config_path.exists():
            raise ValueError(f"Configuration file not found: {config_file}")

        logger.info(f"Loading configuration from {config_file}")

        try:
            with open(config_path, "r") as f:
                # Parse TOML and convert to dictionary
                self._toml_doc = toml_parse(f.read())

                # Parse global patterns
                self._global_patterns = Cfg.parse_global_patterns(self._toml_doc)

                # Get command-specific section or empty dict if not found
                # Command name is derived from the class name (e.g., CollapseCommand -> collapse)
                section = self.__class__.__name__.lower().replace("command", "")

                # Get the section data directly from the TOML document
                if section in self._toml_doc:
                    # Convert the TOML table to a dictionary
                    section_data = {}
                    for key, value in self._toml_doc[section].items():
                        section_data[key] = value
                    self._config = section_data
                    logger.debug(f"Loaded configuration for section '{section}'")
                else:
                    self._config = {}
                    logger.debug(f"No configuration section found for '{section}'")

                return self._config

        except Exception as e:
            logger.error(f"Failed to load configuration from {config_file}: {e}")
            raise ValueError(f"Invalid configuration file {config_file}: {e}")

    def _parse_patterns_from_config(self) -> List[Pattern]:
        """Parse Pattern objects from the loaded configuration.

        Returns
        -------
        List[Pattern]
            List of Pattern objects parsed from configuration

        Raises
        ------
        ValueError
            If no configuration is loaded or pattern configuration is invalid
        """
        if not self._config:
            raise ValueError("No configuration loaded. Call _load_config() first.")

        try:
            patterns = Cfg.parse_patterns_from_config(
                self._config, self._global_patterns
            )
            logger.debug(f"Parsed {len(patterns)} patterns from configuration")
            return patterns
        except Exception as e:
            logger.error(f"Failed to parse patterns from configuration: {e}")
            raise

    def _merge_config_and_args(self, defaults: Optional[Dict[str, Any]] = None) -> None:
        """Merge command line arguments, config file settings, and defaults.

        This method implements a priority system for argument resolution:
        1. Command line arguments (if explicitly set) - highest priority
        2. Config file settings - medium priority
        3. Default values - lowest priority

        Parameters
        ----------
        defaults : Optional[Dict[str, Any]], default=None
            Optional dictionary of default values

        Notes
        -----
        The method tracks which arguments were explicitly set via command line
        to avoid overriding them with config file values.
        """
        if not self.args:
            logger.debug("No arguments to merge")
            return

        # Track which args were explicitly set via command line
        # We do this by checking if the value differs from the default
        if defaults:
            for key, default_value in defaults.items():
                if hasattr(self.args, key):
                    current_value = getattr(self.args, key)
                    if current_value != default_value:
                        self._explicit_args.add(key)
                        logger.debug(
                            f"Argument '{key}' explicitly set via command line"
                        )

        # Start with defaults
        merged = defaults.copy() if defaults else {}

        # First apply defaults to any unset values
        for key, value in merged.items():
            if not hasattr(self.args, key):
                setattr(self.args, key, value)
                logger.debug(f"Applied default value for '{key}': {value}")

        # Then apply config values, but only if the value wasn't explicitly set via command line
        if self._config:
            for key, value in self._config.items():
                if key not in self._explicit_args:
                    setattr(self.args, key, value)
                    logger.debug(f"Applied config value for '{key}': {value}")

    def _progbar_iterable(self, iterable: Iterable, **kwargs) -> Iterable:
        """Wrap an iterable with a progress bar if requested.

        Parameters
        ----------  
        iterable : Iterable
            Iterable to wrap with a progress bar
        kwargs : dict
            Keyword arguments to pass to tqdm

        Returns
        -------
        Iterable
            Iterable wrapped with a progress bar if requested
        """
        if self.args.max_reads is not None:
            iterable = islice(iterable, self.args.max_reads)

        # If progress bar is requested, wrap the iterable with a progress bar
        if self.args.progress_bar:
            return tqdm(iterable, **kwargs)
        else:
            return iterable

    def _chk_exists(self, filenames: Union[str, List[str]]) -> None:
        """Check that files exist and are valid.

        This method validates that all specified files exist and are actual files
        (not directories). It provides clear error messages and exits with
        appropriate error codes if validation fails.

        Parameters
        ----------
        filenames : Union[str, List[str]]
            Single filename or list of filenames to validate

        Raises
        ------
        RuntimeError
            If any filename is a directory
        SystemExit
            If any files don't exist (exits with code 1)
        """
        # Convert single filename to list for uniform processing
        if isinstance(filenames, str):
            filenames = [filenames]

        not_exists = []

        for name in filenames:
            if not os.path.exists(name):
                not_exists.append(name)
                logger.error(f"File does not exist: {name}")
            elif os.path.isdir(name):
                logger.error(f"Path is a directory, not a file: {name}")
                raise RuntimeError(f"Directory not supported: {name}")
            else:
                assert os.path.isfile(name)
                logger.debug(f"File validated: {name}")

        if not_exists:
            logger.error(f"Validation failed: {len(not_exists)} file(s) not found")
            for name in not_exists:
                print(f"DOES NOT EXIST: {name}")
            print(f"EXITING: {len(not_exists)} FILE(S) NOT EXIST")
            sys_exit(1)


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
