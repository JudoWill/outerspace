#!/usr/bin/env python3
"""testing config"""

from outerspace.config import Cfg
from tests.pkgtest.utils import get_filename


def test_config_from_compiled_config():
    """testing config"""
    filename = 'tests/configs/compiled_config.toml'
    cfg = Cfg(filename)
    doc = cfg.read_file()
    
    assert 'read_regxlist' in doc['findseq']
    assert 'read1_regxlist' in doc['findseq']
    assert 'read2_regxlist' in doc['findseq']

    assert 'columns' in doc['collapse']
    assert 'mismatches' in doc['collapse']
    assert 'method' in doc['collapse']

    assert 'barcode_column' in doc['count']
    assert 'key_column' in doc['count']

def test_config_from_grnaquery():
    """testing config"""
    filename = 'tests/configs/grnaquery.toml'
    cfg = Cfg(filename)
    doc = cfg.read_file()
    
    assert 'read_regxlist' in doc['findseq']
    assert doc['findseq']['read1_regxlist'] == [
    "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}",
    "(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?P<downstreamof_protospacer>GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}",
    ]
    
    assert 'read2_regxlist' in doc['findseq']
    assert doc['findseq']['read2_regxlist'] == [
    '(?P<UMI_3prime>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}',
]

    assert 'columns' in doc['collapse']
    assert doc['collapse']['columns'] == 'UMI_5prime,UMI_3prime'
    assert 'mismatches' in doc['collapse']
    assert doc['collapse']['mismatches'] == 2
    assert 'method' in doc['collapse']
    assert doc['collapse']['method'] == 'directional'

    assert 'barcode_column' in doc['count']
    assert doc['count']['barcode_column'] == 'UMI_5prime_UMI_3prime_corrected'
    assert 'key_column' in doc['count']
    assert doc['count']['key_column'] == 'protospacer'

    assert 'column' in doc['gini']
    assert doc['gini']['column'] == 'count'

def test_config_has_not_changed():
    """testing config"""
    filename = 'tests/configs/compiled_config.toml'
    cfg = Cfg(filename)
    doc = cfg.read_file()
    
    default_doc = Cfg.get_doc_default()
    differences = []
    
    # Compare each section
    for section in doc:
        if section not in default_doc:
            differences.append(f"Section '{section}' exists in config but not in default")
            continue
            
        # Compare each field in the section
        for field in doc[section]:
            if field not in default_doc[section]:
                differences.append(f"Field '{field}' in section '{section}' exists in config but not in default")
            elif doc[section][field] != default_doc[section][field]:
                differences.append(f"Field '{field}' in section '{section}' has different value: {doc[section][field]} != {default_doc[section][field]}")
                
    # Check for sections in default that aren't in config
    for section in default_doc:
        if section not in doc:
            differences.append(f"Section '{section}' exists in default but not in config")
            
    assert not differences, "Config differences found:\n" + "\n".join(differences)
    
