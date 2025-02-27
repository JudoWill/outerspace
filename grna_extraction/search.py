"""Search for motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from itertools import chain
from collections import namedtuple
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
        # Getting list and compiling it into object that contains various pieces including  compiled state machine - this is expensive so we are doing it once and then saving
        self.cmps = [regex_compile(p) for p in self.regxlist]
        # Storing names - multi dimensional list- input
        # TODO: Add functionality & tests for multiply nested capture patterns
        self.names = self.capture_names()
        # This is flattened into a list - output - names & patterns
        self.nto = namedtuple('GroupNames', list(chain.from_iterable(self.names)))

    #this is doing the search on the reads using our compiled (sequence) patterns
    def get_capture_from_read(self, read):
        """Searching for provided patterns in the read sequence"""
        ret = []
        # seq is from biopython not a variable
        # assert len(names) == len(self.cmps)
        print('')
        sequence = str(read.seq)
        # 2 groups of names and 2 compiled objects being zipped together
        # Makig sure lists are the same length so that if length are unequal the longer one would not be dropped
        assert len(self.names) == len(self.cmps)
        for idx, (name, cmp) in enumerate(zip(self.names, self.cmps)):
            # TODO: Revisit user interface bc works for us but not world ?Why are we only returning one sequence if look for all sequences?
            # Findall if there were less or more than 1 match we didn't want them- expect only one

            result = cmp.findall(sequence)
            print(f'FINDALL  RESULT{idx}: {result}')
            if len(result) == 1:
                ret.append(self._name_to_captured(name, result[0]))
                print(f'RETURN   RESULT{idx}:  {result[0]}')
        return ret


    # creating empty dictionary
    # adding key values to dictionary name : result
    # zipping together the name and result if result has 2 values or more; Ex. protospacer 1 & 2
    # TODO: what if capture patterns are nested
    def _name_to_captured(self, name, result):
        ret = {}
        print(f'MATCHING:{name}, RESULT: {result}')
        if isinstance(result, str):
            assert len(name) == 1
            ret[name[0]] = result
        else:
            for nam, res in zip(name, result):
                ret[nam] = res
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
            # print(f'{idx}EHHHH, {cmp.groupindex} {cmp.named_lists} {cmp.pattern}')
            # looping cmp.groupindex into name so researcher provided capture group names then get added
            #TODO: want to use their index for list
            names.append(list(cmp.groupindex.keys()))
        return names


    # def _init_nto(self):



# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
