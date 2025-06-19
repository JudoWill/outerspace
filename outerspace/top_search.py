"""Manages multiple search objects for paired reads and non-paired reads"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from outerspace.search import Search


# creating class that manages/holds search objects
class TopSearch:
    """Manages multiple search objects for paired reads and non-paired reads"""

    # srch holds patterns for single read or both paired reads, here searching in both regxlist and list one or list 2
    # srch1 hold patterns for read 1 or paired read
    # read hold patterns for read 2 or paired read
    def __init__(self, read_regxlist=None, read1_regxlist=None, read2_regxlist=None):
        # TODO: Include regxlist in searches for read1 & read2
        self.read_srch = Search(read_regxlist)
        self.read1_srch = Search(read1_regxlist)
        self.read2_srch = Search(read2_regxlist)

    def get_capture_from_readpair(self, read1, read2):
        """Obtaining captured pattern from both reads"""
        match1 = self.read1_srch.get_capture_from_read(read1)
        # TODO: Allow for case where match1 or match2 do not have patterns
        # If match1 exists, adding read id & return match1 
        if match1:
            # If match1 & match2 exists, add match2 to match1
            match2 = self.read2_srch.get_capture_from_read(read2)
            if match2:
                for key,val in match2.items():
                    match1[key] = val
                self._check_readids(read1.id, read2.id)
                match1['read_id'] = read1.id
                return match1
        # TODO: We are only dong match all right now
        # TODO: Add match any
        return None
    #TODO: Think about worlds case when they want partial matches returned, ex. want match any or match all: Rachel=match all

    def get_capture_from_singleread(self, read):
        """Obtaining captured pattern from both reads"""
        mtch = self.read_srch.get_capture_from_read(read)
        if mtch:
            mtch['read_id'] = read.id
        return mtch

    def get_names_readpair(self):
        """Get capture names from search objects for read1 & read2"""
        # TODO: consider self.srch_names
        return ['read_id',] + self.read1_srch.names + self.read2_srch.names  

    def get_names_singleread(self):
        """Get capture names from read search object"""
        # TODO: consider self.srch_names
        return ['read_id',] + self.read_srch.names

    @staticmethod
    def _check_readids(read1_id, read2_id):
        assert read1_id == read2_id, f'ERROR: read1 id not equal to read2 id \n read1: {read1_id}\n read2: {read2_id}'
        return True

    # read_desc is actual name given
    def get_objsearch(self, read_desc):
        """Get requested search object for paired reads or single reads"""
        if read_desc == 'read1':
            # return object if have names if not don't
            return self.read1_srch if self.read1_srch.names else None

        if read_desc == 'read2':
            # return object if have names if not don't
            return self.read if self.read.names else None

        if read_desc == 'read':
            # return object if have names if not don't
            return self.read_srch if self.read_srch.names else None

        # Letting us know if there is something we are trying to run that has not been implemented
        raise RuntimeError(f'UNRECOGNIZED INPUT:{read_desc}: Expected read1, read2, or read')


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
