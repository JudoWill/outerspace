"""Base class for all CLI commands"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
from argparse import ArgumentParser, Namespace
from sys import exit as sys_exit
from typing import Dict, Any, Optional, Set, Type, List
from tomlkit import parse as toml_parse
from outerspace.config import Cfg
from outerspace.pattern import Pattern

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
            self._toml_doc = toml_parse(f.read())
            
            # Parse global patterns
            self._global_patterns = Cfg.parse_global_patterns(self._toml_doc)
            
            # Get command-specific section or empty dict if not found
            # Command name is derived from the class name (e.g., CollapseCommand -> collapse)
            section = self.__class__.__name__.lower().replace('command', '')
            
            # Get the section data directly from the TOML document
            if section in self._toml_doc:
                # Convert the TOML table to a dictionary
                section_data = {}
                for key, value in self._toml_doc[section].items():
                    section_data[key] = value
                self._config = section_data
            else:
                self._config = {}
            
            return self._config

    def _parse_patterns_from_config(self) -> List[Pattern]:
        """Parse Pattern objects from the loaded config.
        
        Returns:
            List of Pattern objects
            
        Raises:
            ValueError: If no config is loaded or pattern configuration is invalid
        """
        if not self._config:
            raise ValueError("No configuration loaded. Call _load_config() first.")
        
        return Cfg.parse_patterns_from_config(self._config, self._global_patterns)

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