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
    assert sorted(dct['define_motifs'].keys()) == sorted(['regex_flags', 'umi_pattern_forward', 'protospacer_pattern', 'back_umi_forward', 'umi_pattern_reverse'])

    
    print(cfg)
    print(dct)
    print(dct['define_motifs'].keys())
    print("TEST PASSED")

if __name__ == '__main__':
    test_config()
