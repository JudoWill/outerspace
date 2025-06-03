"""Main CLI entry point for outerspace"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from argparse import ArgumentParser, Namespace
from sys import exit as sys_exit

from outerspace.cli.commands.findseq import FindSeqCommand
from outerspace.cli.commands.collapse import CollapseCommand
from outerspace.cli.commands.count import CountCommand
from outerspace.cli.commands.gini import GiniCommand
from outerspace.cli.commands.visualize import VisualizeCommand
from outerspace.cli.commands.pipeline import PipelineCommand

class Cli:
    """The command line interface argument requirements"""
    def __init__(self, args=None):
        self.parser = self._init_parser()
        if isinstance(args, Namespace):
            self.args = args
        else:
            self.args = self.parser.parse_args(args)
        self.command = self._init_command()

    def _init_parser(self):
        """Initialize the main argument parser with subcommands"""
        parser = ArgumentParser(
            prog='outerspace',
            description='OUTERSPACE - Optimized Utilities for Tracking Enrichment in Screens',
            epilog='Created by ThreeBlindMice - See how they run code'
        )
        subparsers = parser.add_subparsers(dest='command', help='Available commands')
        
        # Register each command
        for cmd_cls in [FindSeqCommand, CollapseCommand, CountCommand, 
                       GiniCommand, VisualizeCommand, PipelineCommand]:
            cmd = cmd_cls()
            cmd._init_parser(subparsers)
        
        return parser

    def _init_command(self):
        """Initialize the selected command"""
        if not self.args.command:
            return None
            
        command_map = {
            'findseq': FindSeqCommand,
            'collapse': CollapseCommand,
            'count': CountCommand,
            'gini': GiniCommand,
            'visualize': VisualizeCommand,
            'pipeline': PipelineCommand
        }
        return command_map[self.args.command](self.args)

    def run(self):
        """Run the selected command"""
        if not self.command:
            self.parser.print_help()
            return
        self.command.run()

def main(args=None):
    """Main entry point for the CLI"""
    cli = Cli(args)
    cli.run() 