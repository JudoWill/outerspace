"""Logging configuration for OUTERSPACE CLI.

This module provides centralized logging configuration for the OUTERSPACE command-line
interface. It supports both console and file logging with customizable levels and
formats.
"""

import logging
import sys
from pathlib import Path
from typing import Optional

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def setup_logging(
    level: int = logging.INFO, log_file: Optional[str] = None
) -> logging.Logger:
    """Set up logging configuration for OUTERSPACE CLI.

    This function configures the root logger with appropriate handlers and formatters.
    It supports both console output (to stderr) and optional file logging.

    Parameters
    ----------
    level : int, default=logging.INFO
        Logging level for all handlers. Common values:
        - logging.DEBUG: Detailed information for debugging
        - logging.INFO: General information about program execution
        - logging.WARNING: Warning messages for potentially problematic situations
        - logging.ERROR: Error messages for serious problems
        - logging.CRITICAL: Critical errors that may prevent the program from running
    log_file : Optional[str], default=None
        Optional path to log file. If provided, logs will be written to both
        console and file. If None, logs only to console.

    Returns
    -------
    logging.Logger
        Configured logger instance for the calling module

    Notes
    -----
    This function uses force=True in basicConfig to override any existing
    logging configuration. It creates the log file directory if it doesn't exist.

    Examples
    --------
    Basic console logging:
        logger = setup_logging()

    Console and file logging:
        logger = setup_logging(level=logging.DEBUG, log_file="outerspace.log")

    Debug level with custom log file:
        logger = setup_logging(level=logging.DEBUG, log_file="logs/debug.log")
    """
    # Create formatter with timestamp, logger name, level, and message
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Create handlers list
    handlers = []

    # Console handler (stderr for better error handling)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(formatter)
    handlers.append(console_handler)

    # File handler if specified
    if log_file:
        log_path = Path(log_file)
        # Create parent directories if they don't exist
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)

    # Configure root logger with all handlers
    logging.basicConfig(
        level=level,
        handlers=handlers,
        force=True,  # Override any existing configuration
    )

    # Create logger for this module
    logger = logging.getLogger(__name__)
    logger.debug("Logging configuration initialized")

    # Log configuration details
    if log_file:
        logger.info(f"Logging to console and file: {log_file}")
    else:
        logger.info("Logging to console only")

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for a specific module or component.

    This is a convenience function to get a properly configured logger
    for any module in the OUTERSPACE package.

    Parameters
    ----------
    name : str
        Logger name, typically __name__ from the calling module

    Returns
    -------
    logging.Logger
        Configured logger instance

    Examples
    --------
    In a module:
        logger = get_logger(__name__)
        logger.info("Module initialized")
    """
    return logging.getLogger(name)


def set_log_level(level: int) -> None:
    """Set the logging level for the root logger.

    This function allows dynamic adjustment of the logging level
    after initial configuration.

    Parameters
    ----------
    level : int
        New logging level to set

    Examples
    --------
    Enable debug logging:
        set_log_level(logging.DEBUG)

    Reduce verbosity:
        set_log_level(logging.WARNING)
    """
    logging.getLogger().setLevel(level)
    logger = logging.getLogger(__name__)
    logger.debug(f"Log level changed to {logging.getLevelName(level)}")


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
