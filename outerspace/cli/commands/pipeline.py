"""Pipeline command for running complete OUTERSPACE workflows.

This module provides the PipelineCommand class for executing complete OUTERSPACE
pipelines using Snakemake. It handles configuration loading, argument parsing,
and workflow execution with comprehensive error handling.
"""

import logging
import os
import shlex
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from argparse import ArgumentParser

import pulp
import snakemake
import yaml

from outerspace.cli.commands.base import BaseCommand

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

# Monkey patch pulp.list_solvers if it doesn't exist
if not hasattr(pulp, "list_solvers"):
    pulp.list_solvers = pulp.listSolvers


class PipelineCommand(BaseCommand):
    """Command for running the complete OUTERSPACE pipeline using Snakemake.

    This command orchestrates the execution of complete OUTERSPACE workflows
    using Snakemake as the workflow engine. It handles configuration management,
    argument parsing, and workflow execution with comprehensive error handling.
    """

    def __init__(self, args: Optional[Any] = None) -> None:
        """Initialize the pipeline command.

        Parameters
        ----------
        args : Optional[Any], default=None
            Parsed command-line arguments
        """
        super().__init__(args=args)

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "pipeline", help="Run the complete OUTERSPACE pipeline using Snakemake"
        )

        # Required arguments
        parser.add_argument(
            "config_file", help="TOML configuration file with search patterns"
        )
        parser.add_argument(
            "snakemake_config", help="YAML configuration file for Snakemake workflow"
        )

        # Optional arguments
        parser.add_argument(
            "--snakemake-args",
            help='Additional arguments to pass to Snakemake (e.g. --snakemake-args="--dry-run --cores 4")',
        )

    def _load_snakemake_config(self, config_file: str) -> Dict[str, Any]:
        """Load Snakemake configuration from YAML file.

        Parameters
        ----------
        config_file : str
            Path to Snakemake configuration file

        Returns
        -------
        Dict[str, Any]
            Loaded configuration dictionary

        Raises
        ------
        ValueError
            If configuration file doesn't exist or is invalid
        """
        logger.info(f"Loading Snakemake config from: {config_file}")

        if not os.path.exists(config_file):
            raise ValueError(f"Snakemake configuration file not found: {config_file}")

        try:
            with open(config_file, "r") as f:
                config = yaml.safe_load(f)
                logger.debug(f"Loaded config: {config}")
                return config
        except Exception as e:
            logger.error(f"Failed to load Snakemake config from {config_file}: {e}")
            raise ValueError(f"Invalid Snakemake configuration file {config_file}: {e}")

    def _parse_snakemake_args(self, args_str: Optional[str]) -> List[str]:
        """Parse additional Snakemake arguments.

        Parameters
        ----------
        args_str : Optional[str]
            String of arguments to parse (e.g. "--dry-run --cores 4")

        Returns
        -------
        List[str]
            List of parsed arguments

        Notes
        -----
        Uses shlex.split to handle quoted arguments correctly and preserve
        argument boundaries.
        """
        if not args_str:
            return []

        # Parse the arguments using shlex to handle quotes correctly
        parsed_args = shlex.split(args_str)
        logger.debug(f"Parsed Snakemake arguments: {parsed_args}")
        return parsed_args

    def run(self) -> None:
        """Run the pipeline command using Snakemake.

        This method orchestrates the pipeline execution by loading configurations,
        preparing Snakemake arguments, and executing the workflow with comprehensive
        error handling and logging.

        Raises
        ------
        ValueError
            If required configuration files are missing or invalid
        SystemExit
            If Snakemake execution fails (exit code != 0)
        """
        logger.info(f"Running pipeline with config file: {self.args.config_file}")
        logger.info(f"Snakemake config file: {self.args.snakemake_config}")

        # Load TOML config
        self._load_config(self.args.config_file)

        # Load Snakemake config
        snakemake_config = self._load_snakemake_config(self.args.snakemake_config)

        # Prepare Snakemake configuration
        config = {"toml": self.args.config_file, **snakemake_config}

        # Parse additional Snakemake arguments
        snakemake_args = self._parse_snakemake_args(self.args.snakemake_args)

        # Construct new argv for Snakemake
        # Start with the program name
        new_argv = ["snakemake"]

        # Add config file if specified
        if "--configfile" not in snakemake_args:
            new_argv.extend(["--configfile", self.args.snakemake_config])

        # Add remaining arguments
        new_argv.extend(snakemake_args)

        # Add our config as a config value
        if "--config" not in snakemake_args:
            new_argv.extend(["--config", f"toml={self.args.config_file}"])

        logger.info(f"Running Snakemake with args: {new_argv}")

        # Run Snakemake using its main function
        try:
            logger.info("Starting Snakemake workflow execution")
            snakemake.main(new_argv[1:])  # Skip the 'snakemake' program name
            logger.info("Pipeline completed successfully")

        except SystemExit as e:
            if e.code != 0:
                logger.error(f"Pipeline failed with exit code {e.code}")
                sys.exit(1)
            else:
                logger.info("Pipeline completed successfully")
        except Exception as e:
            logger.error(f"Unexpected error during pipeline execution: {e}")
            raise


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
