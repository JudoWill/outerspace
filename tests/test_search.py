#!/usr/bin/env python3
"""testing config"""

from grna_extraction.config import Cfg
from grna_extraction.search import Search
print(f'TTTTTTTT test_search({__name__})')


def test_search():
    """testing config"""
    cfg = Cfg()
    doc = cfg.get_doc_default()
    regxlist = doc['regxlist']
    search = Search(regxlist)
    print(regxlist)
    print("TEST PASSED")

if __name__ == '__main__':
    test_search()
