"""Pipeline command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
from pathlib import Path
import sys
from typing import Dict, Any, List
import yaml
import logging
import shlex

from outerspace.cli.commands.base import BaseCommand
import snakemake

import pulp

# Monkey patch pulp.list_solvers if it doesn't exist
if not hasattr(pulp, 'list_solvers'):
    pulp.list_solvers = pulp.listSolvers

logger = logging.getLogger(__name__)

class PipelineCommand(BaseCommand):
    """Command for running the complete OUTERSPACE pipeline using Snakemake"""
    
    def __init__(self, args=None):
        super().__init__(args=args)
    
    def _init_parser(self, subparsers):
        """Initialize command-specific argument parser"""
        parser = subparsers.add_parser('pipeline',
            help='Run the complete OUTERSPACE pipeline using Snakemake')
        
        # Required arguments
        parser.add_argument('config_file',
            help='TOML configuration file with search patterns')
        parser.add_argument('snakemake_config',
            help='YAML configuration file for Snakemake workflow')
        
        # Optional arguments
        parser.add_argument('--snakemake-args',
            help='Additional arguments to pass to Snakemake (e.g. --snakemake-args="--dry-run --cores 4")')
        
        return parser

    def _load_snakemake_config(self, config_file: str) -> Dict[str, Any]:
        """Load Snakemake configuration from YAML file"""
        logger.info(f"Loading Snakemake config from: {config_file}")
        if not os.path.exists(config_file):
            raise ValueError(f"Snakemake configuration file not found: {config_file}")
        
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
            logger.debug(f"Loaded config: {config}")
            return config

    def _parse_snakemake_args(self, args_str: str) -> List[str]:
        """Parse additional Snakemake arguments
        
        Args:
            args_str: String of arguments to parse (e.g. "--dry-run --cores 4")
            
        Returns:
            List of parsed arguments
        """
        if not args_str:
            return []
            
        # Parse the arguments using shlex to handle quotes correctly
        return shlex.split(args_str)

    def run(self):
        """Run the pipeline command using Snakemake"""
        logger.info(f"Running pipeline with config file: {self.args.config_file}")
        logger.info(f"Snakemake config file: {self.args.snakemake_config}")
        
        # Load TOML config
        self._load_config(self.args.config_file)
        
        # Load Snakemake config
        snakemake_config = self._load_snakemake_config(self.args.snakemake_config)
        
        # Prepare Snakemake configuration
        config = {
            'toml': self.args.config_file,
            **snakemake_config
        }
        
        # Parse additional Snakemake arguments
        snakemake_args = self._parse_snakemake_args(self.args.snakemake_args)
        
        # Construct new argv for Snakemake
        # Start with the program name
        new_argv = ['snakemake']
        
        # Add config file if specified
        if '--configfile' not in snakemake_args:
            new_argv.extend(['--configfile', self.args.snakemake_config])
        
        # Add remaining arguments
        new_argv.extend(snakemake_args)
        
        # Add our config as a config value
        if '--config' not in snakemake_args:
            new_argv.extend(['--config', f"toml={self.args.config_file}"])
        
        logger.info(f"Running Snakemake with args: {new_argv}")
        
        # Run Snakemake using its main function
        try:
            snakemake.main(new_argv[1:])  # Skip the 'snakemake' program name
        except SystemExit as e:
            if e.code != 0:
                logger.error("Pipeline failed")
                sys.exit(1)
        
        logger.info("Pipeline completed successfully") 