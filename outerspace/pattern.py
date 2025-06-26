"""Pattern matching functionality for sequence analysis.

This module provides classes for defining and executing search patterns on DNA/RNA
sequences. It includes the Pattern class for pattern definition and the Hit class
for storing match results.
"""
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import logging
from typing import Dict, Generator, List, Optional, Union
from regex import compile as regex_compile, Pattern as RegexPattern

from outerspace.read import Read

# Set up logging
logger = logging.getLogger(__name__)


class Hit:
    """Container for pattern match information.

    This class stores the details of a successful pattern match, including
    position information, the matched sequence, orientation, and any captured
    groups from the regular expression.
    """

    def __init__(
        self,
        start: int,
        end: int,
        match: str,
        orientation: str,
        captured: Dict[str, str],
    ) -> None:
        """Initialize a Hit object with match information.

        Parameters
        ----------
        start : int
            Starting position of the match in the sequence
        end : int
            Ending position of the match in the sequence (exclusive)
        match : str
            The actual matched sequence
        orientation : str
            Orientation of the match ('forward' or 'reverse-complement')
        captured : Dict[str, str]
            Dictionary of named capture groups from the regex
        """
        self.start = start
        self.end = end
        self.match = match
        self.orientation = orientation
        self.captured = captured

    def __str__(self) -> str:
        """Return string representation of the hit.

        Returns
        -------
        str
            Formatted string showing hit details
        """
        return (
            f"Hit(start={self.start}, end={self.end}, match={self.match}, "
            f"orientation={self.orientation}, captured={self.captured})"
        )

    def __repr__(self) -> str:
        """Return detailed string representation for debugging.

        Returns
        -------
        str
            Same as __str__ for this class
        """
        return self.__str__()


class Pattern:
    """Pattern definition for sequence searching.

    This class defines a search pattern with a regular expression and search
    parameters including read selection, orientation, and multiple match handling.
    """

    # Valid values for pattern parameters
    VALID_READS = {"R1", "R2", "both"}
    VALID_ORIENTATIONS = {"forward", "reverse-complement", "both"}
    VALID_MULTIPLE_OPTIONS = {"first", "last", "all"}

    def __init__(
        self,
        reg_expr: str,
        read: str,
        orientation: str,
        multiple: str,
        search_read_name: bool = False,
        left_flank: int = 0,
        right_flank: int = 0,
    ) -> None:
        """Initialize a Pattern object.

        Parameters
        ----------
        reg_expr : str
            Regular expression pattern to search for
        read : str
            Which read to search ('R1', 'R2', or 'both')
        orientation : str
            Search orientation ('forward', 'reverse-complement', or 'both')
        multiple : str
            How to handle multiple matches ('first', 'last', or 'all')
        search_read_name : bool
            Whether to search the read name for the pattern
        left_flank : int
            Number of bases to include before the match
        right_flank : int
            Number of bases to include after the match

        Raises
        ------
        ValueError
            If any parameter has an invalid value
        """
        self.reg_expr = reg_expr
        self.read = read
        self.orientation = orientation
        self.multiple = multiple
        self.search_read_name = search_read_name
        self.left_flank = left_flank
        self.right_flank = right_flank

        # Validate parameters
        self._check_args()

        # Compile the regular expression
        self._regex = self._regex_compile()

        logger.debug(
            f"Initialized pattern: {reg_expr} on {read} read, "
            f"{orientation} orientation, {multiple} matches"
        )

    def _check_args(self) -> None:
        """Validate pattern parameters.

        Raises
        ------
        ValueError
            If any parameter has an invalid value
        """
        if self.read not in self.VALID_READS:
            raise ValueError(
                f"Invalid read: {self.read}. Must be one of {self.VALID_READS}"
            )

        if self.orientation not in self.VALID_ORIENTATIONS:
            raise ValueError(
                f"Invalid orientation: {self.orientation}. "
                f"Must be one of {self.VALID_ORIENTATIONS}"
            )

        if self.multiple not in self.VALID_MULTIPLE_OPTIONS:
            raise ValueError(
                f"Invalid multiple: {self.multiple}. "
                f"Must be one of {self.VALID_MULTIPLE_OPTIONS}"
            )

        if not isinstance(self.search_read_name, bool):
            raise ValueError(
                f"Invalid search_read_name: {self.search_read_name}. "
                f"Must be a boolean"
            )

        if not isinstance(self.left_flank, int) or not isinstance(
            self.right_flank, int
        ):
            raise ValueError(
                f"Invalid flanks: {self.left_flank} and {self.right_flank}. "
                f"Must be integers"
            )

    def __str__(self) -> str:
        """Return string representation of the pattern.

        Returns
        -------
        str
            Formatted string showing pattern parameters
        """
        return (
            f"Pattern(reg_expr={self.reg_expr}, read={self.read}, "
            f"orientation={self.orientation}, multiple={self.multiple}, "
            f"search_read_name={self.search_read_name}, "
            f"left_flank={self.left_flank}, right_flank={self.right_flank})"
        )

    def __repr__(self) -> str:
        """Return detailed string representation for debugging.

        Returns
        -------
        str
            Same as __str__ for this class
        """
        return self.__str__()

    def _regex_compile(self) -> RegexPattern:
        """Compile the regular expression pattern.

        Returns
        -------
        RegexPattern
            Compiled regular expression object

        Notes
        -----
        Uses the 'regex' library instead of 're' for enhanced functionality.
        """
        try:
            return regex_compile(self.reg_expr)
        except Exception as e:
            logger.error(f"Failed to compile regex pattern '{self.reg_expr}': {e}")
            raise

    @staticmethod
    def _search(
        regex: RegexPattern,
        sequence: str,
        orientation: str,
        left_flank: int,
        right_flank: int,
    ) -> Generator[Hit, None, None]:
        """Search for matches in a sequence.

        Parameters
        ----------
        regex : RegexPattern
            Compiled regular expression to search with
        sequence : str
            Sequence to search in
        orientation : str
            Orientation label for the matches
        left_flank : int
            Number of bases to include before the match
        right_flank : int
            Number of bases to include after the match

        Yields
        ------
        Hit
            Hit objects for each match found
        """
        for match in regex.finditer(sequence):
            start = max(0, match.start() - left_flank)
            end = min(len(sequence), match.end() + right_flank)
            match_str = sequence[start:end]
            yield Hit(
                start=start,
                end=end,
                match=match_str,
                orientation=orientation,
                captured=match.groupdict(),
            )

    def _search_read(self, read: Read) -> Generator[Hit, None, None]:
        """Search for pattern matches in a read.

        This method handles read selection and orientation searching based on
        the pattern configuration.

        Parameters
        ----------
        read : Read
            Read object to search in

        Yields
        ------
        Hit
            Hit objects for each match found
        """

        # Search the read name for the pattern
        if self.search_read_name:
            yield from self._search(
                self._regex, read.name, "forward", self.left_flank, self.right_flank
            )
        else:
            # Check if this read should be searched
            if (self.read == "both") or (read.pair == self.read):
                # Search forward orientation
                if (self.orientation == "forward") | (self.orientation == "both"):
                    yield from self._search(
                        self._regex,
                        read.seq,
                        "forward",
                        self.left_flank,
                        self.right_flank,
                    )

                # Search reverse-complement orientation
                if (self.orientation == "reverse-complement") | (
                    self.orientation == "both"
                ):
                    yield from self._search(
                        self._regex,
                        read.seq_rc,
                        "reverse-complement",
                        self.left_flank,
                        self.right_flank,
                    )

    def search(self, read: Read) -> Optional[Union[Hit, List[Hit]]]:
        """Search for pattern matches in a read.

        This method returns matches based on the 'multiple' configuration:
        - 'first': Returns the first match found or None
        - 'last': Returns the last match found or None
        - 'all': Returns a list of all matches (empty list if none found)

        Parameters
        ----------
        read : Read
            Read object to search in

        Returns
        -------
        Optional[Union[Hit, List[Hit]]]
            Match results based on the 'multiple' configuration:
            - Single Hit object for 'first' or 'last'
            - List of Hit objects for 'all'
            - None if no matches found for 'first' or 'last'
            - Empty list if no matches found for 'all'
        """
        try:
            if self.multiple == "first":
                result = next(self._search_read(read))
                logger.debug(f"Found first match: {result}")
                return result
            elif self.multiple == "last":
                matches = list(self._search_read(read))
                if matches:
                    result = matches[-1]
                    logger.debug(f"Found last match: {result}")
                    return result
                else:
                    logger.debug("No matches found for 'last' search")
                    return None
            else:  # 'all'
                matches = list(self._search_read(read))
                logger.debug(f"Found {len(matches)} matches for 'all' search")
                return matches

        except StopIteration:
            logger.debug("No matches found for 'first' search")
            return None
        except Exception as e:
            logger.error(f"Error during pattern search: {e}")
            raise


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
