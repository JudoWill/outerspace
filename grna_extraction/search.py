"""Search for motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from regex import compile as regex_compile


# creating class that will do the searching
class Search:
    """Search for motifs"""

    def __init__(self, regxlist=None):
        """making the regxlist accessible so that all the functions inside this class can see it"""
        self.regxlist = regxlist if regxlist is not None else []
        # p is pattern
        # comprehension utilized below instead of for loop across multiple lines
        # evertime loop through p will add an element 
        # cmp= compile
        # TODO: ?dictionary of compile object with names?
        self.cmps = [regex_compile(p) for p in self.regxlist]
    
    #this is doing the search on the reads using our compiled (sequence) patterns
    def get_capture_from_read(self, read):
        """Searching for provided patterns in the read sequence"""
        # seq is from biopython not a variable     
        ret = []
        sequence = str(read.seq)
        for cmp in self.cmps:
            # TODO: ?Why are we only returning one sequence if look for all sequences?
            result = cmp.findall(sequence)
            if len(result) == 1:
                ret.append(result[0])
        return ret
    


    # Want to write a function that captures names from regxlist
    # Name is within ex. <UMI> 
    # https://pypi.org/project/regex/
    # names of functions/datamembers returned included: 'findall', 'finditer', 'flags', 'fullmatch', 'groupindex', 'groups', 'match', 'named_lists', 'pattern', 'scanner'
    def capture_names(self):
        """Get names from regxlist to use as column headers"""
        # IMPORTANT FOR USER EXPERIENCE:
        # TODO: Account for 2 capture names being the same in one regex pattern
        # TODO: Account for 2 capture names being the same across regex pattern
        # names is empty list, if dictionary would be {}
        names = []     
        # for index & compile object enumerate the objects
        for idx, cmp in enumerate(self.cmps):
            #print(dir(cmp))
            # search empty string ""
            #mtch = cmp.search("")
            print(f'{idx}EHHHH, {cmp.groupindex} {cmp.named_lists} {cmp.pattern}')
            for e in cmp.groupindex:
                print(f'{idx}AHHH({e})')
        return names




# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
