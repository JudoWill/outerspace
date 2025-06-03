"""Base class for all CLI commands"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
from argparse import ArgumentParser, Namespace
from sys import exit as sys_exit
from typing import Dict, Any, Optional, Set
from tomlkit import parse as toml_parse

class BaseCommand:
    """Base class for all commands"""
    def __init__(self, args=None):
        self.args = args
        self._config = None
        self._explicit_args: Set[str] = set()

    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        raise NotImplementedError("Each command must implement _init_parser")

    def run(self):
        """Execute the command"""
        raise NotImplementedError("Each command must implement run")

    def _load_config(self, config_file: str) -> Dict[str, Any]:
        """Load settings from TOML configuration file"""
        if not os.path.exists(config_file):
            raise ValueError(f"Configuration file not found: {config_file}")
        
        with open(config_file, 'r') as f:
            # Parse TOML and convert to dictionary
            toml_doc = toml_parse(f.read())
            
            # Get command-specific section or empty dict if not found
            # Command name is derived from the class name (e.g., CollapseCommand -> collapse)
            section = self.__class__.__name__.lower().replace('command', '')
            
            # Get the section data directly from the TOML document
            if section in toml_doc:
                # Convert the TOML table to a dictionary
                section_data = {}
                for key, value in toml_doc[section].items():
                    section_data[key] = value
                self._config = section_data
            else:
                self._config = {}
            
            return self._config

    def _merge_config_and_args(self, defaults: Optional[Dict[str, Any]] = None) -> None:
        """Merge command line arguments, config file settings, and defaults into self.args.
        
        Priority order (highest to lowest):
        1. Command line arguments (if explicitly set)
        2. Config file settings
        3. Default values
        
        Args:
            defaults: Optional dictionary of default values
        """
        if not self.args:
            return

        # Track which args were explicitly set via command line
        # We do this by checking if the value differs from the default
        if defaults:
            for key, default_value in defaults.items():
                if hasattr(self.args, key):
                    current_value = getattr(self.args, key)
                    if current_value != default_value:
                        self._explicit_args.add(key)

        # Start with defaults
        merged = defaults.copy() if defaults else {}
        
        # First apply defaults to any unset values
        for key, value in merged.items():
            if not hasattr(self.args, key):
                setattr(self.args, key, value)
        
        # Then apply config values, but only if the value wasn't explicitly set via command line
        if self._config:
            for key, value in self._config.items():
                if key not in self._explicit_args:
                    setattr(self.args, key, value)

    def _chk_exists(self, filenames):
        """Check that files exist"""
        from os.path import exists, isdir, isfile
        not_exists = []
        for name in filenames:
            if not exists(name):
                not_exists.append(name)
            elif isdir(name):
                raise RuntimeError(f'Directory not supported: {name}')
            else:
                assert isfile(name)
        if not_exists:
            for name in not_exists:
                print(f'DOES NOT EXIST: {name}')
            print(f'EXITING: {len(not_exists)} FILE(S) NOT EXIST')
            sys_exit(1) 