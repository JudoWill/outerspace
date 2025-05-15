#!/usr/bin/env python3
"""testing config"""

from outerspace.config import Cfg
from tests.pkgtest.utils import get_filename
print(f'TTTTTTTT test_config({__name__})')


def test_config():
    """testing config"""
    filename = get_filename('example.cfg')
    cfg = Cfg(filename)
    doc = cfg.write_file()
    dct = doc.unwrap()
    assert doc is not None 
    # assert not dct, f"Dictionary contains stuff & shouldn't:\n{dct}"
    assert dct, f"Dictionary is empty & should contain stuff: {dct}"
    assert sorted(dct.keys()) == ['owner', 'regxlist', 'regxlist1', 'regxlist2', 'title'], sorted(dct.keys())

    return
    # TODO: Adjust tests below
    # Asserting variable is equal to what is defined in outerspace/config.py
    assert cfg.get_umi_pattern_forward() == (
        '(?P<UMI>.{8})'
        '(?:CTTGGCTTTATATATCTTGTGG)'
        '{s<=4}')
    print(cfg.get_umi_pattern_forward())

    # Asserting variable is equal to what is defined in outerspace/config.py
    assert cfg.get_protospacer_forward() == (
        '(?:TATCTTGTGGAAAGGACGAAACACC)'
        '{s<=4}'
        '(?P<protospacer>.{19,21})'
        '(?:GTTTAAGTACTCTGTGCTGGAAACAG)'
        '{s<=4}')
    print(cfg.get_protospacer_forward())

    #Asserting variable is equal to what is defined in outerspace/config.py
    exp_pat = 'gtgtgtcagttagggtgtggaa'
    assert cfg.get_umi_pattern_forward_downstream_nt() == exp_pat,(
            f'\nACTUAL: {cfg.get_umi_pattern_forward_downstream_nt()}\n'
            f'EXPECTED: {exp_pat}')
            
    print(cfg.get_umi_pattern_forward_downstream_nt())

    #Asserting variable is equal to what is defined in outerspace/config.py
    

        
    # These prints below were to see the keys that were printed out to then add them to the assert tests above
    # print(cfg)
    # print(dct)
    # print(dct['define_motifs'].keys())

    print("TEST PASSED")

if __name__ == '__main__':
    test_config()
