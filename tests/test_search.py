#!/usr/bin/env python3
"""testing config"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"

from grna_extraction.config import Cfg
from grna_extraction.search import Search
print(f'TTTTTTTT test_search({__name__})')


def test_search():
    """testing config"""
    cfg = Cfg()
    doc = cfg.get_doc_default()
    regxlist = doc['regxlist']
    search = Search(regxlist)
    names = search.capture_names()
    print(regxlist)
    print(names)
    print("TEST PASSED")
    

if __name__ == '__main__':
    test_search()





 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
