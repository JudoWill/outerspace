#!/usr/bin/env python3
"""testing config"""

import pytest
from outerspace.config import Cfg
from tests.pkgtest.utils import get_filename
import tempfile
import os
from outerspace.pattern import Pattern

@pytest.mark.xfail(reason="Config has changed. #TODO update config")
def test_config_from_compiled_config():
    """testing config"""
    filename = "tests/configs/compiled_config.toml"
    cfg = Cfg(filename)
    doc = cfg.read_file()

    # Check for new global patterns structure
    assert "patterns" in doc
    assert len(doc["patterns"]) > 0

    # Check that patterns have names
    for pattern in doc["patterns"]:
        assert "name" in pattern

    assert "columns" in doc["collapse"]
    assert "mismatches" in doc["collapse"]
    assert "method" in doc["collapse"]

    assert "barcode_column" in doc["count"]
    assert "key_column" in doc["count"]


def test_config_from_grnaquery():
    """testing config"""
    filename = "tests/configs/grnaquery.toml"
    cfg = Cfg(filename)
    doc = cfg.read_file()

    # Check for new global patterns structure
    assert "patterns" in doc
    assert len(doc["patterns"]) == 3

    # Check pattern names
    pattern_names = [p["name"] for p in doc["patterns"]]
    assert "UMI_5prime" in pattern_names
    assert "protospacer" in pattern_names
    assert "UMI_3prime" in pattern_names

    # Check findseq section uses pattern names
    assert "pattern_names" in doc["findseq"]
    assert doc["findseq"]["pattern_names"] == [
        "UMI_5prime",
        "protospacer",
        "UMI_3prime",
    ]

    assert "columns" in doc["collapse"]
    assert doc["collapse"]["columns"] == "UMI_5prime,UMI_3prime"
    assert "mismatches" in doc["collapse"]
    assert doc["collapse"]["mismatches"] == 2
    assert "method" in doc["collapse"]
    assert doc["collapse"]["method"] == "directional"

    assert "barcode_column" in doc["count"]
    assert doc["count"]["barcode_column"] == "UMI_5prime_UMI_3prime_corrected"
    assert "key_column" in doc["count"]
    assert doc["count"]["key_column"] == "protospacer"


@pytest.mark.xfail(reason="Config has changed")
def test_config_has_not_changed():
    """testing config"""
    filename = "tests/configs/compiled_config.toml"
    cfg = Cfg(filename)
    doc = cfg.read_file()

    default_doc = Cfg.get_doc_default()
    differences = []

    # Compare each section
    for section in doc:
        if section not in default_doc:
            differences.append(
                f"Section '{section}' exists in config but not in default"
            )
            continue

        # Compare each field in the section
        for field in doc[section]:
            if field not in default_doc[section]:
                differences.append(
                    f"Field '{field}' in section '{section}' exists in config but not in default"
                )
            elif doc[section][field] != default_doc[section][field]:
                differences.append(
                    f"Field '{field}' in section '{section}' has different value: {doc[section][field]} != {default_doc[section][field]}"
                )

    # Check for sections in default that aren't in config
    for section in default_doc:
        if section not in doc:
            differences.append(
                f"Section '{section}' exists in default but not in config"
            )

    assert not differences, "Config differences found:\n" + "\n".join(differences)


def test_parse_new_pattern_format():
    """Test parsing new pattern format from config"""
    config_data = {
        "patterns": [
            {
                "name": "UMI_5prime",
                "reg_expr": "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}",
                "read": "R1",
                "orientation": "forward",
                "multiple": "first",
            },
            {
                "name": "protospacer",
                "reg_expr": "(?P<protospacer>.{19,21})",
                "read": "R1",
                "orientation": "both",
                "multiple": "all",
            },
            {
                "name": "UMI_3prime",
                "reg_expr": "(?P<UMI_3prime>.{8})",
                "read": "R2",
                "orientation": "reverse-complement",
                "multiple": "last",
            },
        ]
    }

    patterns = Cfg.parse_patterns_from_config(config_data)

    assert len(patterns) == 3

    # Check first pattern
    assert (
        patterns[0].reg_expr == "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}"
    )
    assert patterns[0].name == "UMI_5prime"
    assert patterns[0].read == "R1"
    assert patterns[0].orientation == "forward"
    assert patterns[0].multiple == "first"

    # Check second pattern
    assert patterns[1].reg_expr == "(?P<protospacer>.{19,21})"
    assert patterns[1].name == "protospacer"
    assert patterns[1].read == "R1"
    assert patterns[1].orientation == "both"
    assert patterns[1].multiple == "all"

    # Check third pattern
    assert patterns[2].reg_expr == "(?P<UMI_3prime>.{8})"
    assert patterns[2].name == "UMI_3prime"
    assert patterns[2].read == "R2"
    assert patterns[2].orientation == "reverse-complement"
    assert patterns[2].multiple == "last"


def test_parse_empty_config():
    """Test parsing empty config"""
    config_data = {}
    patterns = Cfg.parse_patterns_from_config(config_data)
    assert len(patterns) == 0


def test_parse_invalid_pattern_values():
    """Test parsing patterns with invalid values"""
    config_data = {
        "patterns": [
            {
                "name": "test",
                "reg_expr": "(?P<test>.{5})",
                "read": "INVALID",  # Invalid read value
                "orientation": "forward",
                "multiple": "first",
            }
        ]
    }

    with pytest.raises(ValueError, match="Invalid read"):
        Cfg.parse_patterns_from_config(config_data)


def test_parse_toml_file_with_patterns():
    """Test parsing patterns from actual TOML file"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".toml", delete=False) as f:
        f.write(
            """[findseq]
[[findseq.patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}"
read = "R1"
orientation = "forward"
multiple = "first"

[[findseq.patterns]]
name = "protospacer"
reg_expr = "(?P<protospacer>.{19,21})"
read = "R1"
orientation = "both"
multiple = "all"
"""
        )
        temp_file = f.name

    try:
        # Read the TOML file
        cfg = Cfg(temp_file)
        config_data = cfg.read_file()

        # Parse patterns
        patterns = Cfg.parse_patterns_from_config(config_data["findseq"])

        assert len(patterns) == 2
        assert (
            patterns[0].reg_expr
            == "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}"
        )
        assert patterns[0].name == "UMI_5prime"
        assert patterns[0].read == "R1"
        assert patterns[1].reg_expr == "(?P<protospacer>.{19,21})"
        assert patterns[1].name == "protospacer"
        assert patterns[1].orientation == "both"

    finally:
        os.unlink(temp_file)


def test_parse_global_patterns():
    """Test parsing global patterns from TOML document"""
    from tomlkit import parse

    toml_str = """# Global patterns
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})"
read = "R1"
orientation = "forward"
multiple = "first"

[[patterns]]
name = "protospacer"
reg_expr = "(?P<protospacer>.{19,21})"
read = "R1"
orientation = "forward"
multiple = "first"
"""

    toml_doc = parse(toml_str)
    global_patterns = Cfg.parse_global_patterns(toml_doc)

    assert len(global_patterns) == 2
    assert "UMI_5prime" in global_patterns
    assert "protospacer" in global_patterns

    # Check pattern content
    umi_pattern = global_patterns["UMI_5prime"]
    assert umi_pattern["reg_expr"] == "(?P<UMI_5prime>.{8})"
    assert umi_pattern["read"] == "R1"


def test_parse_patterns_with_global_patterns():
    """Test parsing patterns using global pattern references"""
    from tomlkit import parse

    toml_str = """# Global patterns
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})"
read = "R1"
orientation = "forward"
multiple = "first"

[[patterns]]
name = "protospacer"
reg_expr = "(?P<protospacer>.{19,21})"
read = "R1"
orientation = "forward"
multiple = "first"

[findseq]
pattern_names = ["UMI_5prime", "protospacer"]
"""

    toml_doc = parse(toml_str)
    global_patterns = Cfg.parse_global_patterns(toml_doc)

    # Parse findseq section
    findseq_config = toml_doc["findseq"]
    patterns = Cfg.parse_patterns_from_config(findseq_config, global_patterns)

    assert len(patterns) == 2
    assert patterns[0].reg_expr == "(?P<UMI_5prime>.{8})"
    assert patterns[1].reg_expr == "(?P<protospacer>.{19,21})"


def test_parse_patterns_with_use_all_patterns():
    """Test parsing patterns using use_all_patterns flag"""
    from tomlkit import parse

    toml_str = """# Global patterns
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})"
read = "R1"
orientation = "forward"
multiple = "first"

[[patterns]]
name = "protospacer"
reg_expr = "(?P<protospacer>.{19,21})"
read = "R1"
orientation = "forward"
multiple = "first"

[findseq]
use_all_patterns = true
"""

    toml_doc = parse(toml_str)
    global_patterns = Cfg.parse_global_patterns(toml_doc)

    # Parse findseq section
    findseq_config = toml_doc["findseq"]
    patterns = Cfg.parse_patterns_from_config(findseq_config, global_patterns)

    assert len(patterns) == 2
    # Order might vary, so check both patterns exist
    reg_exprs = [p.reg_expr for p in patterns]
    assert "(?P<UMI_5prime>.{8})" in reg_exprs
    assert "(?P<protospacer>.{19,21})" in reg_exprs


def test_parse_patterns_missing_global_pattern():
    """Test error handling for missing global pattern"""
    from tomlkit import parse

    toml_str = """# Global patterns
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})"
read = "R1"
orientation = "forward"
multiple = "first"

[findseq]
pattern_names = ["UMI_5prime", "missing_pattern"]
"""

    toml_doc = parse(toml_str)
    global_patterns = Cfg.parse_global_patterns(toml_doc)

    # Parse findseq section
    findseq_config = toml_doc["findseq"]

    with pytest.raises(
        ValueError, match="Pattern 'missing_pattern' not found in global patterns"
    ):
        Cfg.parse_patterns_from_config(findseq_config, global_patterns)


def test_parse_global_patterns_missing_name():
    """Test error handling for global patterns without name"""
    from tomlkit import parse

    toml_str = """# Global patterns
[[patterns]]
reg_expr = "(?P<UMI_5prime>.{8})"
read = "R1"
orientation = "forward"
multiple = "first"
"""

    toml_doc = parse(toml_str)

    with pytest.raises(ValueError, match="Global patterns must have a 'name' field"):
        Cfg.parse_global_patterns(toml_doc)
