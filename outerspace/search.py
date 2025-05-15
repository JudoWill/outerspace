"""Search for motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "SC Barrera"

from itertools import chain
from collections import namedtuple
from logging import debug
from regex import compile as regex_compile



# creating class that will do the searching
class Search:
    """Search for motifs"""

    def __init__(self, regxlist=None):
        """making the regxlist accessible so that all the functions inside this class can see it"""
        # saving a copy of the orignial list of strings
        self.regxlist = regxlist if regxlist is not None else []
        # p is pattern
        # comprehension utilized below instead of for loop across multiple lines
        # evertime loop through p will add an element
        # cmp= compile
        # TODO: ?dictionary of compile object with names?
        # Creating a list and compiling it into object that contains various pieces including  compiled state machine - this is expensive so we are doing it once and then saving
        self.cmps = [regex_compile(p) for p in self.regxlist]
        # Storing names - multi dimensional list- input
        # TODO: Add functionality & tests for multiply nested capture patterns
        self.capturednames = self._init_capturednames()
        # This is flattened into a list - output - names & patterns
        # self.nto = namedtuple('GroupNames', list(chain.from_iterable(self.capturednames)))
        # Make a set once of the names from each compile bc change from iterable is expensive
        # this is flat list that retains order
        self.names = list(chain.from_iterable(self.capturednames))
        # This is set for for checking existance bc sets do not retain order
        self.nameset = set(self.names)


    def get_csvheaders(self):
        """These are the headers for the csv containing output search results"""
        # This is what is being written into csv
        return ['read_id'] + self.names

    #this is doing the search on the reads using our compiled (sequence) patterns
    def get_capture_from_read(self, read):
        """Searching for provided patterns in the read sequence"""
        # if no capture patterns found don't return an empty dictionary
        if not self.names:
            return None
        #TODO FOR RACHEL: add read2
        read_dict = {}
        # seq is from biopython not a variable
        # assert len(names) == len(self.cmps)
        # debug is used instead of print but can turn it on or off from any file, ex turned off in test
        debug('')
        sequence = str(read.seq)
        # 2 groups of names and 2 compiled objects being zipped together
        # Makig sure lists are the same length so that if length are unequal the longer one would not be dropped
        assert len(self.capturednames) == len(self.cmps)
        for idx, (name, cmp) in enumerate(zip(self.capturednames, self.cmps)):
            # TODO FIRST: Add findall option to check we don't have multiple matches, search gives back first match
            # TODO: Revisit user interface bc works for us but not world ?Why are we only returning one sequence if look for all sequences?
            # Search if there was a pattern or not- looking for the first pattern
            #Match utiliing search
            mtch = cmp.search(sequence)
            # print(f'SEARCH  RESULT{idx}: {mtch}')
            # If search for pattern successful then return object; if no pattern found none is returned 
            if mtch is not None:
                #adding search items found (key & value) into dictionary
                for key, val in mtch.groupdict().items():
                    read_dict[key] = val
                # read_dict.append(self._name_to_captured(name, mtch)
                # print(f'SEARCHDICT      cmp{idx}:  {mtch.groupdict()}')
                # print(f'SEARCH   GROUP0-cmp{idx}:  {mtch.group(0)}')
                # print(f'SEARCH   GROUP1-cmp{idx}:  {mtch.group(1)}')
        debug(f'CHECKED DICTIONARY:  {read_dict}')        
        # TODO:we are only reporting if we found everything but that is not what the world will expect- we need to report if anything was found 
        # set is used to check between 2 lists; printing if something was returned
        if self.nameset == set(read_dict.keys()):
            # read_dict['read_id'] = read.id
            debug(f'RETURN DICTIONARY: {read_dict}')
            return read_dict
        return None

    # This is for findall not search
    # Findall if there were less or more than 1 match we didn't want them- expect only one
    def _run_findall_(self, ret, idx, cmp, sequence, name):
        # TODO FIRST: Add findall option to check we don't have multiple matches, search gives back first match
        find_list = cmp.findall(sequence)
        print(f'FINDALL  RESULT{idx}: {find_list}')
        if len(find_list) == 1:
            ret.append(self._name_to_captured(name, find_list[0]))
            print(f'RETURN   RESULT{idx}:  {find_list[0]}')
        return ret


    # creating empty dictionary
    # adding key values to dictionary name : result
    # zipping together the name and result if result has 2 values or more; Ex. protospacer 1 & 2
    # TODO: what if capture patterns are nested
    @staticmethod
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
    def _init_capturednames(self):
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
