"""Hold search patterns """

from regex import compile as regex_compile

from outerspace.read import Read

class Hit:
    """Hold hit information"""

    def __init__(self, start: int, end: int, match: str, orientation: str, captured: dict):
        self.start = start
        self.end = end
        self.match = match
        self.orientation = orientation
        self.captured = captured

    def __str__(self):
        return f"Hit(start={self.start}, end={self.end}, match={self.match}, orientation={self.orientation}, captured={self.captured})"

    def __repr__(self):
        return self.__str__()


class Pattern:
    """Hold search patterns"""

    def __init__(self, reg_expr: str, read: str, orientation: str, multiple: str):
        self.reg_expr = reg_expr
        self.read = read
        self.orientation = orientation
        self.multiple = multiple
        self._check_args()

        self._regex = self._regex_compile()
        

    def _check_args(self):
        if self.read not in ['R1', 'R2', 'both']:
            raise ValueError(f"Invalid read: {self.read}")
        if self.orientation not in ['forward', 'reverse-complement', 'both']:
            raise ValueError(f"Invalid orientation: {self.orientation}")
        if self.multiple not in ['first', 'last', 'all']:
            raise ValueError(f"Invalid multiple: {self.multiple}")

    def __str__(self):
        return f"Pattern(reg_expr={self.reg_expr}, read={self.read}, orientation={self.orientation}, multiple={self.multiple})"

    def __repr__(self):
        return self.__str__()

    def _regex_compile(self):
        return regex_compile(self.reg_expr)

    @staticmethod
    def _search(regex, sequence: str, orientation: str):
        for match in regex.finditer(sequence):
            yield Hit(
                start=match.start(),
                end=match.end(),
                match=match.group(0),
                orientation=orientation,
                captured=match.groupdict()
            )
    
    def _search_read(self, read: Read):
        """Return all hits"""
        if (self.read == 'both') or (read.pair == self.read):
            # Search forward orientation
            if self.orientation in ['forward', 'both']:
                yield from self._search(self._regex, read.seq, 'forward')
            
            # Search reverse-complement orientation
            if self.orientation in ['reverse-complement', 'both']:
                yield from self._search(self._regex, read.seq_rc, 'reverse-complement')
    
    def search(self, read):
        """Return first, last, or all hits"""
        
        try:
            if self.multiple == 'first':
                return next(self._search_read(read))
            elif self.multiple == 'last':
                return list(self._search_read(read))[-1]
            else:
                return list(self._search_read(read))
        except StopIteration:
            return None
