"""Manages multiple search objects for paired reads and non-paired reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from grna_extraction.search import Search


# creating class that manages/holds search objects
class TopSearch:
    """Manages multiple search objects for paired reads and non-paired reads"""

    # srch holds patterns for single read or both paired reads
    # srch1 hold patterns for read 1 or paired read
    # srch2 hold patterns for read 2 or paired read
    def __init__(self, regxdct=None):
        self.srch = Search(regxdct['regxlist'])
        self.srch1 = Search(regxdct['regxlist1'])
        self.srch2 = Search(regxdct['regxlist2'])

    # read_desc is actual name given
    def get_objsearch(self, read_desc):
        """Get requested search object for paired reads or single reads"""
        #TODO: Get search for single read or both paired reads- not yet needed for Rachel
        if read_desc == 'read1':
            # assert there is nothing in the list, print out what is in the list if there is
            assert not self.srch.cmps, self.srch.regxlist
            return self.srch1
        if read_desc == 'read2':
            # assert there is nothing in the list, print out what is in the list if there is
            assert not self.srch.cmps, self.srch.regxlist
            return self.srch2
        # Letting us know if there is something we are trying to run that has not been implemented
        raise NotImplementedError(f'UNRECOGNIZED INPUT: TopSearch.get_objsearch({read_desc}')
    

# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
