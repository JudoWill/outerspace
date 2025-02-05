#!/usr/bin/env python3
"""testing config"""

from grna_extraction.config import Cfg
print(f'TTTTTTTT test_config({__name__})')


def test_config():
    """testing config"""
    cfg = Cfg()
    dct = cfg.doc.unwrap()
    assert cfg.doc is not None 
    # assert not dct, f"Dictionary contains stuff & shouldn't:\n{dct}"
    assert dct, f"Dictionary is empty & should contain stuff: {dct}"
    assert sorted(dct.keys()) == sorted(['title', 'owner', 'define_motifs'])
    assert sorted(dct['define_motifs'].keys()) == sorted([
        'regex_flags',
        "umi_pattern_forward_num_umi_nt",
        "umi_pattern_forward_pattern_nt",
        "umi_pattern_forward_mismatch_max_nt",
        "pattern_forward_upstream_protospacer",
        "pattern_forward_upstream_protospacer_mismatch_max_nt",
        "pattern_forward_num_protospacer_nt_range_from",
        "pattern_forward_num_protospacer_nt_range_to",
        "pattern_forward_downstream_protospacer",
        "pattern_forward_protospacer_mismatch_max_nt",
        "umi_pattern_forward_downstream_nt",
        'umi_pattern_reverse'])

    # Asserting variable is equal to what is defined in grna_extraction/config.py
    assert cfg.get_umi_pattern_forward() == (
        '(?P<UMI>.{8})'
        '(?:CTTGGCTTTATATATCTTGTGG)'
        '{s<=4}')
    print(cfg.get_umi_pattern_forward())

    # Asserting variable is equal to what is defined in grna_extraction/config.py
    assert cfg.get_protospacer_forward() == (
        '(?:TATCTTGTGGAAAGGACGAAACACC)'
        '{s<=4}'
        '(?P<protospacer>.{19,21})'
        '(?:GTTTAAGTACTCTGTGCTGGAAACAG)'
        '{s<=4}')
    print(cfg.get_protospacer_forward())

    #Asserting variable is equal to what is defined in grna_extraction/config.py
    assert cfg.get_umi_pattern_forward_downstream_nt() == ('gtgtgtcagttagggtgtggaa')
    print(cfg.get_umi_pattern_forward_downstream_nt())

    #Asserting variable is equal to what is defined in grna_extraction/config.py
    

        
    # These prints below were to see the keys that were printed out to then add them to the assert tests above
    # print(cfg)
    # print(dct)
    # print(dct['define_motifs'].keys())

    print("TEST PASSED")

if __name__ == '__main__':
    test_config()
