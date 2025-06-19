"""Logging configuration for OUTERSPACE CLI"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import logging
import sys
from pathlib import Path

def setup_logging(level=logging.INFO, log_file=None):
    """Set up logging configuration
    
    Args:
        level: Logging level (default: INFO)
        log_file: Optional path to log file (default: None, logs to stderr)
    """
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Create handlers
    handlers = []
    
    # Console handler (stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(formatter)
    handlers.append(console_handler)
    
    # File handler if specified
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        handlers.append(file_handler)
    
    # Configure root logger
    logging.basicConfig(
        level=level,
        handlers=handlers,
        force=True  # Override any existing configuration
    )
    
    # Create logger for this module
    logger = logging.getLogger(__name__)
    logger.debug("Logging configured")
    
    return logger 