"""Manages multiple search objects for paired reads and non-paired reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from grna_extraction.search import Search


# creating class that manages/holds search objects
class TopSearch:
    """Manages multiple search objects for paired reads and non-paired reads"""

    # srch holds patterns for single read or both paired reads, here searching in both regxlist and list one or list 2
    # srch1 hold patterns for read 1 or paired read
    # srch2 hold patterns for read 2 or paired read
    def __init__(self, regxdct=None):
        self.srch = Search(regxdct['regxlist'])
        self.srch1 = Search(regxdct['regxlist'] + regxdct['regxlist1'])
        self.srch2 = Search(regxdct['regxlist'] + regxdct['regxlist2'])

    def get_capture_from_readpair(self, read1, read2):
        """Obtaining captured pattern from both reads"""
        match1 = self.srch1.get_capture_from_read(read1)
        match2 = self.srch2.get_capture_from_read(read2)
        # If math1 exists and match2 exists, FININSH THIS
        if match1:
            if match2:
                for key,val in match2.items():
                    match1[key] = val
            return match1
        return match2

    def get_capture_from_singleread(self, read):
        """Obtaining captured pattern from both reads"""
        return self.srch.get_capture_from_read(read)

    # read_desc is actual name given
    def get_objsearch(self, read_desc):
        """Get requested search object for paired reads or single reads"""
        if read_desc == 'read1':
            # return object if have names if not don't
            return self.srch1 if self.srch1.names else None

        if read_desc == 'read2':
            # return object if have names if not don't
            return self.srch2 if self.srch2.names else None

        if read_desc == 'read':
            # return object if have names if not don't
            return self.srch if self.srch.names else None

        # Letting us know if there is something we are trying to run that has not been implemented
        raise RuntimeError(f'UNRECOGNIZED INPUT:{read_desc}: Expected read1, read2, or read')


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
