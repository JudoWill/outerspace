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
        # TODO: Include regxlist in searches for read1 & read2
        self.srch = Search(regxdct['regxlist'])
        self.srch1 = Search(regxdct['regxlist1'])
        self.srch2 = Search(regxdct['regxlist2'])
        print(f'captured names for read1: {self.srch1.capturednames}')
        print(f'captured names for read2: {self.srch2.capturednames}')

    def get_capture_from_readpair(self, read1, read2):
        """Obtaining captured pattern from both reads"""
        match1 = self.srch1.get_capture_from_read(read1)
        # TODO: Allow for case where match1 or match2 do not have patterns
        # If match1 exists, adding read id & return match1 
        if match1:
            # If match1 & match2 exists, add match2 to match1
            match2 = self.srch2.get_capture_from_read(read2)
            if match2:
                for key,val in match2.items():
                    match1[key] = val
                self._check_readids(read1.id, read2.id)
                match1['read_id'] = read1.id
                return match1
        # TODO: We are only dong match all right now
        # TODO: Add match any
        # if match1 does not exist does match2 exist, if so return match2
        # if match2:
           # match2['read_id'] = read2.id
        # return match2
        return None
    #TODO: Think about worlds case when they want partial matches returned, ex. want match any or match all: Rachel=match all

    def get_capture_from_singleread(self, read):
        """Obtaining captured pattern from both reads"""
        mtch = self.srch.get_capture_from_read(read)
        if mtch:
            mtch['read_id'] = read.id
        return mtch

    def get_names_readpair(self):
        """Get capture names from search objects for read1 & read2"""
        # TODO: consider self.srch_names
        return ['read_id',] + self.srch1.names + self.srch2.names  

    def get_names_singleread(self):
        """Get capture names from read search object"""
        # TODO: consider self.srch_names
        return ['read_id',] + self.srch.names

    @staticmethod
    def _check_readids(read1_id, read2_id):
        if read1_id == read2_id:
            return True
        # TODO: Write this to a log file
        # TODO: Summarize mismatched reads to STDOUT
        print(f'ERROR: read1 id not equal to read2 id \n read1: {read1_id}\n read2: {read2_id}')
        return False

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
