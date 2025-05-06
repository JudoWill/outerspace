#!/usr/bin/env python3
"""testing config"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"

from grna_extraction.config import Cfg
from grna_extraction.search import Search
print(f'TTTTTTTT test_search({__name__})')


def test_search():
    """testing config"""
    # Creating a config object; allows researcher to provide search patterns
    # Using default configuration file
    # TODO: may want a default configuration in examples
    # 

    cfg = Cfg()
    doc = cfg.get_doc_default()
    regxlist = doc['regxlist1']
    search = Search(regxlist)
    names = search.capturednames
    assert regxlist
    print(f'CAPTURED REGXLIST CONSISTS OF: {regxlist}')
    assert names, f'EXPECTED NAMES, GOT: {names}'
    assert names == [['UMI'], ['protospacer', 'downstreamof_protospacer']], names
    print(f' NAMES: {search.names}')
   
    print("TEST PASSED")
    

if __name__ == '__main__':
    test_search()





 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
