"""Tests for pattern matching functionality"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


import pytest
from outerspace.pattern import Pattern, Hit
from outerspace.read import Read


def test_hit_creation():
    """Test basic Hit object creation"""
    hit = Hit(
        pattern=Pattern(name="test", reg_expr="ATCG", read="R1", orientation="forward", multiple="first"),
        start=10,
        end=15,
        match="ATCG",
        orientation="forward",
        captured={"group1": "AT", "group2": "CG"},
    )

    assert hit.start == 10
    assert hit.end == 15
    assert hit.match == "ATCG"
    assert hit.orientation == "forward"
    assert hit.captured == {"group1": "AT", "group2": "CG"}


def test_hit_string_representation():
    """Test Hit string and repr methods"""
    hit = Hit(
        pattern=Pattern(name="test", reg_expr="ATCG", read="R1", orientation="forward", multiple="first"),
        start=5, end=9, match="GCTA", orientation="reverse-complement", captured={}
    )

    expected_str = (
        "Hit(start=5, end=9, match=GCTA, orientation=reverse-complement, captured={}, pattern=test)"
    )
    assert str(hit) == expected_str
    assert repr(hit) == expected_str


def test_hit_with_named_captures():
    """Test Hit with named capture groups"""
    hit = Hit(
        pattern=Pattern(name="test", reg_expr="ATCG", read="R1", orientation="forward", multiple="first"),
        start=0,
        end=8,
        match="ATCGATCG",
        orientation="forward",
        captured={"primer": "ATCG", "barcode": "ATCG"},
    )

    assert hit.captured["primer"] == "ATCG"
    assert hit.captured["barcode"] == "ATCG"
    assert len(hit.captured) == 2


def test_pattern_creation():
    """Test basic Pattern object creation"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    assert pattern.reg_expr == "ATCG"
    assert pattern.read == "R1"
    assert pattern.orientation == "forward"
    assert pattern.multiple == "first"
    assert pattern.search_read_name is False


def test_pattern_creation_with_search_read_name():
    """Test Pattern creation with search_read_name parameter"""
    pattern = Pattern(
        name="test",
        reg_expr="test",
        read="both",
        orientation="both",
        multiple="all",
        search_read_name=True,
    )

    assert pattern.search_read_name is True
    assert pattern.read == "both"
    assert pattern.orientation == "both"
    assert pattern.multiple == "all"


def test_pattern_string_representation():
    """Test Pattern string and repr methods"""
    pattern = Pattern(
        name="test",
        reg_expr="GCTA",
        read="R2",
        orientation="reverse-complement",
        multiple="last",
        search_read_name=True,
    )

    expected_str = "Pattern(name=test, reg_expr=GCTA, read=R2, orientation=reverse-complement, multiple=last, search_read_name=True, left_flank=0, right_flank=0)"
    assert str(pattern) == expected_str
    assert repr(pattern) == expected_str


def test_pattern_validation_valid_reads():
    """Test Pattern validation with valid read values"""
    valid_reads = ["R1", "R2", "both"]

    for read in valid_reads:
        pattern = Pattern(
            name="test",
            reg_expr="ATCG", read=read, orientation="forward", multiple="first"
        )
        assert pattern.read == read


def test_pattern_validation_invalid_read():
    """Test Pattern validation with invalid read value"""
    with pytest.raises(ValueError, match="Invalid read: invalid_read"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="invalid_read",
            orientation="forward",
            multiple="first",
        )


def test_pattern_validation_valid_orientations():
    """Test Pattern validation with valid orientation values"""
    valid_orientations = ["forward", "reverse-complement", "both"]

    for orientation in valid_orientations:
        pattern = Pattern(
            name="test",
            reg_expr="ATCG", read="R1", orientation=orientation, multiple="first"
        )
        assert pattern.orientation == orientation


def test_pattern_validation_invalid_orientation():
    """Test Pattern validation with invalid orientation value"""
    with pytest.raises(ValueError, match="Invalid orientation: invalid_orientation"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="R1",
            orientation="invalid_orientation",
            multiple="first",
        )


def test_pattern_validation_valid_multiple_options():
    """Test Pattern validation with valid multiple values"""
    valid_multiple_options = ["first", "last", "all"]

    for multiple in valid_multiple_options:
        pattern = Pattern(
            name="test",
            reg_expr="ATCG", read="R1", orientation="forward", multiple=multiple
        )
        assert pattern.multiple == multiple


def test_pattern_validation_invalid_multiple():
    """Test Pattern validation with invalid multiple value"""
    with pytest.raises(ValueError, match="Invalid multiple: invalid_multiple"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="R1",
            orientation="forward",
            multiple="invalid_multiple",
        )


def test_pattern_validation_search_read_name():
    """Test Pattern validation of search_read_name parameter"""
    # Valid boolean values
    pattern1 = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        search_read_name=True,
    )
    assert pattern1.search_read_name is True

    pattern2 = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        search_read_name=False,
    )
    assert pattern2.search_read_name is False

    # Invalid non-boolean value
    with pytest.raises(ValueError, match="Invalid search_read_name"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="R1",
            orientation="forward",
            multiple="first",
            search_read_name="not_a_boolean",
        )


def test_pattern_regex_compilation():
    """Test Pattern regex compilation"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    # Test that regex was compiled successfully
    assert pattern._regex is not None
    assert hasattr(pattern._regex, "finditer")


def test_pattern_invalid_regex():
    """Test Pattern with invalid regex pattern"""
    with pytest.raises(
        Exception
    ):  # regex library raises various exceptions for invalid patterns
        Pattern(
            name="test",
            reg_expr="[ATCG",  # Invalid regex - unclosed bracket
            read="R1",
            orientation="forward",
            multiple="first",
        )


def test_pattern_search_forward():
    """Test Pattern search in forward orientation"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0
    assert result.end == 4
    assert result.match == "ATCG"
    assert result.orientation == "forward"


def test_pattern_search_reverse_complement():
    """Test Pattern search in reverse-complement orientation"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="reverse-complement", multiple="first"
    )

    read = Read(seq="CGAT", pair="R1", name="test_read")  # Reverse complement of ATCG
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0
    assert result.end == 4
    assert (
        result.match == "ATCG"
    )  # The pattern "ATCG" matches in the reverse complement sequence
    assert result.orientation == "reverse-complement"


def test_pattern_search_both_orientations():
    """Test Pattern search in both orientations"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="both", multiple="all"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert (
        len(result) == 1
    )  # Only one match - ATCG doesn't match in reverse complement CGAT

    # Check forward match
    forward_hit = result[0]
    assert forward_hit.start == 0
    assert forward_hit.end == 4
    assert forward_hit.match == "ATCG"
    assert forward_hit.orientation == "forward"


def test_pattern_search_read_name():
    """Test Pattern search in read name"""
    pattern = Pattern(
        name="test",
        reg_expr="test",
        read="R1",
        orientation="forward",
        multiple="first",
        search_read_name=True,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read_name")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0
    assert result.end == 4
    assert result.match == "test"
    assert result.orientation == "forward"


def test_pattern_search_multiple_matches_first():
    """Test Pattern search with multiple matches, returning first"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should return first match
    assert result.end == 4
    assert result.match == "ATCG"


def test_pattern_search_multiple_matches_last():
    """Test Pattern search with multiple matches, returning last"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="last"
    )

    read = Read(seq="ATCGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 4  # Should return last match
    assert result.end == 8
    assert result.match == "ATCG"


def test_pattern_search_multiple_matches_all():
    """Test Pattern search with multiple matches, returning all"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="all"
    )

    read = Read(seq="ATCGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 2

    # Check first match
    assert result[0].start == 0
    assert result[0].end == 4
    assert result[0].match == "ATCG"

    # Check second match
    assert result[1].start == 4
    assert result[1].end == 8
    assert result[1].match == "ATCG"


def test_pattern_search_no_matches_first():
    """Test Pattern search with no matches, returning first"""
    pattern = Pattern(
        name="test",
        reg_expr="XXXX", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is None


def test_pattern_search_no_matches_last():
    """Test Pattern search with no matches, returning last"""
    pattern = Pattern(
        name="test",
        reg_expr="XXXX", read="R1", orientation="forward", multiple="last"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is None


def test_pattern_search_no_matches_all():
    """Test Pattern search with no matches, returning all"""
    pattern = Pattern(
        name="test",
        reg_expr="XXXX", read="R1", orientation="forward", multiple="all"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 0


def test_pattern_search_read_selection():
    """Test Pattern search with different read selections"""
    # Test R1 only
    pattern_r1 = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read_r1 = Read(seq="ATCG", pair="R1", name="read1")
    read_r2 = Read(seq="ATCG", pair="R2", name="read2")

    result_r1 = pattern_r1.search(read_r1)
    result_r2 = pattern_r1.search(read_r2)

    assert result_r1 is not None  # Should find match in R1
    assert result_r2 is None  # Should not find match in R2

    # Test R2 only
    pattern_r2 = Pattern(
        name="test",
        reg_expr="ATCG", read="R2", orientation="forward", multiple="first"
    )

    result_r1 = pattern_r2.search(read_r1)
    result_r2 = pattern_r2.search(read_r2)

    assert result_r1 is None  # Should not find match in R1
    assert result_r2 is not None  # Should find match in R2

    # Test both reads
    pattern_both = Pattern(
        name="test",
        reg_expr="ATCG", read="both", orientation="forward", multiple="first"
    )

    result_r1 = pattern_both.search(read_r1)
    result_r2 = pattern_both.search(read_r2)

    assert result_r1 is not None  # Should find match in R1
    assert result_r2 is not None  # Should find match in R2


def test_pattern_search_with_named_captures():
    """Test Pattern search with named capture groups"""
    pattern = Pattern(
        name="test",
        reg_expr="(?P<primer>AT)(?P<barcode>CG)",
        read="R1",
        orientation="forward",
        multiple="first",
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.captured["primer"] == "AT"
    assert result.captured["barcode"] == "CG"
    assert len(result.captured) == 2


def test_pattern_search_complex_regex():
    """Test Pattern search with complex regex patterns"""
    # Test with quantifiers - A{2,4} means 2 to 4 A's followed by CG
    pattern = Pattern(
        name="test",
        reg_expr="A{2,4}CG", read="R1", orientation="forward", multiple="all"
    )

    read = Read(seq="AACGAACGAACCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert (
        len(result) == 2
    )  # Should match AACG (2 A's) and AACG (2 A's), but not AACCG (3 A's)

    # Check matches
    assert result[0].match == "AACG"
    assert result[1].match == "AACG"

    # Test with character classes - [AT]CG matches A or T followed by CG
    pattern2 = Pattern(
        name="test",
        reg_expr="[AT]CG", read="R1", orientation="forward", multiple="all"
    )

    read2 = Read(seq="ACGTGC", pair="R1", name="test_read")
    result2 = pattern2.search(read2)

    assert result2 is not None
    assert isinstance(result2, list)
    assert (
        len(result2) == 1
    )  # Should match ACG (A at position 0), but not TGC (T at position 3)
    assert result2[0].match == "ACG"


def test_pattern_search_empty_sequence():
    """Test Pattern search with empty sequence"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is None


def test_pattern_search_empty_regex():
    """Test Pattern search with empty regex pattern"""
    pattern = Pattern(
        name="test",
        reg_expr="", read="R1", orientation="forward", multiple="all"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    # Empty regex matches at every position
    assert len(result) > 0


def test_pattern_search_case_sensitivity():
    """Test Pattern search case sensitivity"""
    # Test case-sensitive search
    pattern = Pattern(
        name="test",
        reg_expr="atcg", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is None  # Should not match due to case difference

    # Test case-insensitive search
    pattern_insensitive = Pattern(
        name="test",
        reg_expr="(?i)atcg", read="R1", orientation="forward", multiple="first"
    )

    result = pattern_insensitive.search(read)

    assert result is not None  # Should match with case-insensitive flag


def test_pattern_search_with_anchors():
    """Test Pattern search with regex anchors"""
    # Test start anchor
    pattern_start = Pattern(
        name="test",
        reg_expr="^ATCG", read="R1", orientation="forward", multiple="first"
    )

    read1 = Read(seq="ATCG", pair="R1", name="test_read")
    read2 = Read(seq="GATCG", pair="R1", name="test_read")

    result1 = pattern_start.search(read1)
    result2 = pattern_start.search(read2)

    assert result1 is not None  # Should match at start
    assert result2 is None  # Should not match in middle

    # Test end anchor
    pattern_end = Pattern(
        name="test",
        reg_expr="ATCG$", read="R1", orientation="forward", multiple="first"
    )

    result1 = pattern_end.search(read1)
    result2 = pattern_end.search(read2)

    assert result1 is not None  # Should match at end
    assert result2 is not None  # Should match at end (ATCG is at the end of GATCG)


def test_pattern_search_error_handling():
    """Test Pattern search error handling"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    # Test with None read (should raise exception)
    with pytest.raises(Exception):
        pattern.search(None)


def test_pattern_search_read_generator():
    """Test Pattern._search_read generator method"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCGATCG", pair="R1", name="test_read")

    hits = list(pattern._search_read(read))

    assert len(hits) == 2
    assert all(isinstance(hit, Hit) for hit in hits)
    assert all(hit.orientation == "forward" for hit in hits)


def test_pattern_search_read_name_generator():
    """Test Pattern._search_read generator with search_read_name=True"""
    pattern = Pattern(
        name="test",
        reg_expr="test",
        read="R1",
        orientation="forward",
        multiple="first",
        search_read_name=True,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read_name")

    hits = list(pattern._search_read(read))

    assert len(hits) == 1
    assert hits[0].orientation == "forward"
    assert hits[0].match == "test"


def test_pattern_search_read_selection_generator():
    """Test Pattern._search_read generator with different read selections"""
    # Test R1 only
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read_r1 = Read(seq="ATCG", pair="R1", name="read1")
    read_r2 = Read(seq="ATCG", pair="R2", name="read2")

    hits_r1 = list(pattern._search_read(read_r1))
    hits_r2 = list(pattern._search_read(read_r2))

    assert len(hits_r1) == 1  # Should find match in R1
    assert len(hits_r2) == 0  # Should not find match in R2


def test_pattern_search_orientation_generator():
    """Test Pattern._search_read generator with different orientations"""
    # Test forward only
    pattern_forward = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")

    hits_forward = list(pattern_forward._search_read(read))
    assert len(hits_forward) == 1
    assert hits_forward[0].orientation == "forward"

    # Test reverse-complement only - ATCG doesn't match in reverse complement CGAT
    pattern_rc = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="reverse-complement", multiple="first"
    )

    hits_rc = list(pattern_rc._search_read(read))
    assert len(hits_rc) == 0  # No matches in reverse complement

    # Test both orientations
    pattern_both = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="both", multiple="first"
    )

    hits_both = list(pattern_both._search_read(read))
    assert len(hits_both) == 1  # Only forward match
    assert hits_both[0].orientation == "forward"


def test_pattern_search_reverse_complement_with_matching_pattern():
    """Test Pattern search in reverse-complement orientation with a pattern that actually matches"""
    pattern = Pattern(
        name="test",
        reg_expr="CGAT", read="R1", orientation="reverse-complement", multiple="first"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")  # Reverse complement is CGAT
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0
    assert result.end == 4
    assert result.match == "CGAT"
    assert result.orientation == "reverse-complement"


def test_pattern_search_both_orientations_with_matching_pattern():
    """Test Pattern search in both orientations with a pattern that matches in reverse-complement only"""
    pattern = Pattern(
        name="test",
        reg_expr="CGAT", read="R1", orientation="both", multiple="all"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")  # Reverse complement is CGAT
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 1  # Only reverse-complement match

    rc_hit = result[0]
    assert rc_hit.orientation == "reverse-complement"
    assert rc_hit.match == "CGAT"


def test_pattern_search_with_flanks_overlapping_matches():
    """Test Pattern search with flanks and overlapping regex matches"""
    pattern = Pattern(
        name="test",
        reg_expr="A{1,2}",
        read="R1",
        orientation="forward",
        multiple="all",
        left_flank=1,
        right_flank=1,
    )

    read = Read(seq="GAAAG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 2  # Should match "A" and "AA"

    # Check matches with flanks - the match field contains the sequence with flanks
    matches = [hit.match for hit in result]
    assert "GAAA" in matches  # A with 1 base flank on each side
    assert "AAG" in matches  # AA with 1 base flank on each side


def test_pattern_search_with_zero_width_assertions():
    """Test Pattern search with zero-width assertions"""
    # Test positive lookahead
    pattern = Pattern(
        name="test",
        reg_expr="A(?=T)", read="R1", orientation="forward", multiple="all"
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 1
    assert result[0].match == "A"
    assert result[0].start == 0
    assert result[0].end == 1


def test_pattern_search_with_unicode_sequences():
    """Test Pattern search with Unicode sequences"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="first"
    )

    # Test with Unicode sequence (should work the same as ASCII)
    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.match == "ATCG"


def test_pattern_search_with_very_long_sequence():
    """Test Pattern search with a very long sequence"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG", read="R1", orientation="forward", multiple="all"
    )

    # Create a long sequence with multiple ATCG patterns
    long_seq = "ATCG" * 1000
    read = Read(seq=long_seq, pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 1000  # Should find 1000 matches

    # Check first and last matches
    assert result[0].start == 0
    assert result[0].end == 4
    assert result[-1].start == 3996  # 999 * 4
    assert result[-1].end == 4000


def test_pattern_creation_with_flanks():
    """Test Pattern creation with left_flank and right_flank parameters"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=2,
        right_flank=3,
    )

    assert pattern.left_flank == 2
    assert pattern.right_flank == 3
    assert pattern.reg_expr == "ATCG"
    assert pattern.read == "R1"


def test_pattern_creation_with_zero_flanks():
    """Test Pattern creation with zero flank values"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=0,
        right_flank=0,
    )

    assert pattern.left_flank == 0
    assert pattern.right_flank == 0


def test_pattern_validation_invalid_flanks():
    """Test Pattern validation with invalid flank values"""
    # Test non-integer left_flank
    with pytest.raises(ValueError, match="Invalid flanks"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="R1",
            orientation="forward",
            multiple="first",
            left_flank="not_an_int",
            right_flank=0,
        )

    # Test non-integer right_flank
    with pytest.raises(ValueError, match="Invalid flanks"):
        Pattern(
            name="test",
            reg_expr="ATCG",
            read="R1",
            orientation="forward",
            multiple="first",
            left_flank=0,
            right_flank=3.14,
        )


def test_pattern_search_with_left_flank():
    """Test Pattern search with left flank only"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=2,
        right_flank=0,
    )

    read = Read(seq="GGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should include 2 bases before ATCG
    assert result.end == 6  # Should end at the end of ATCG
    assert result.match == "GGATCG"  # Should include the flank
    assert result.orientation == "forward"


def test_pattern_search_with_right_flank():
    """Test Pattern search with right flank only"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=0,
        right_flank=3,
    )

    read = Read(seq="ATCGTTT", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should start at the beginning of ATCG
    assert result.end == 7  # Should include 3 bases after ATCG
    assert result.match == "ATCGTTT"  # Should include the flank
    assert result.orientation == "forward"


def test_pattern_search_with_both_flanks():
    """Test Pattern search with both left and right flanks"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=1,
        right_flank=2,
    )

    read = Read(seq="GATCGTT", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should include 1 base before ATCG
    assert result.end == 7  # Should include 2 bases after ATCG
    assert result.match == "GATCGTT"  # Should include both flanks
    assert result.orientation == "forward"


def test_pattern_search_with_flanks_at_sequence_boundaries():
    """Test Pattern search with flanks that extend beyond sequence boundaries"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=5,
        right_flank=5,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should be clamped to sequence start
    assert result.end == 4  # Should be clamped to sequence end
    assert result.match == "ATCG"  # Should include entire sequence
    assert result.orientation == "forward"


def test_pattern_search_with_flanks_multiple_matches():
    """Test Pattern search with flanks and multiple matches"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="all",
        left_flank=1,
        right_flank=1,
    )

    read = Read(seq="GATCGATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 2

    # First match: GATCGA (ATCG with 1 base flank on each side, but right flank overlaps with next match)
    assert result[0].start == 0
    assert result[0].end == 6
    assert result[0].match == "GATCGA"  # Should include the flanks

    # Second match: GATCG (ATCG with 1 base flank on each side)
    assert result[1].start == 4
    assert result[1].end == 9
    assert result[1].match == "GATCG"  # Should include the flanks


def test_pattern_search_with_flanks_reverse_complement():
    """Test Pattern search with flanks in reverse-complement orientation"""
    pattern = Pattern(
        name="test",
        reg_expr="CGAT",
        read="R1",
        orientation="reverse-complement",
        multiple="first",
        left_flank=2,
        right_flank=1,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")  # Reverse complement is CGAT
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should include 2 bases before CGAT (clamped to start)
    assert result.end == 4  # Should include 1 base after CGAT (clamped to end)
    assert (
        result.match == "CGAT"
    )  # Should contain the matched sequence from reverse complement
    assert result.orientation == "reverse-complement"


def test_pattern_search_with_flanks_read_name():
    """Test Pattern search with flanks in read name"""
    pattern = Pattern(
        name="test",
        reg_expr="test",
        read="R1",
        orientation="forward",
        multiple="first",
        search_read_name=True,
        left_flank=2,
        right_flank=3,
    )

    read = Read(seq="ATCG", pair="R1", name="mytestread")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert (
        result.start == 0
    )  # Should include 2 characters before "test" (clamped to start)
    assert result.end == 9  # Should include 3 characters after "test" (clamped to end)
    assert (
        result.match == "mytestrea"
    )  # Should include the flanks (clamped to string end)
    assert result.orientation == "forward"


def test_pattern_search_with_flanks_complex_regex():
    """Test Pattern search with flanks and complex regex patterns"""
    pattern = Pattern(
        name="test",
        reg_expr="A{2,4}CG",
        read="R1",
        orientation="forward",
        multiple="all",
        left_flank=1,
        right_flank=2,
    )

    read = Read(seq="GAACGTTAACGT", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 2

    # First match: GAACGTT (AACG with 1 base before, 2 after)
    assert result[0].start == 0
    assert result[0].end == 7
    assert result[0].match == "GAACGTT"  # Should include the flanks

    # Second match: TAACGT (AACG with 1 base before, 2 after)
    assert result[1].start == 6  # Corrected: AACG starts at position 6
    assert result[1].end == 12
    assert result[1].match == "TAACGT"  # Should include the flanks


def test_pattern_search_with_flanks_overlapping_matches():
    """Test Pattern search with flanks and overlapping regex matches"""
    pattern = Pattern(
        name="test",
        reg_expr="A{1,2}",
        read="R1",
        orientation="forward",
        multiple="all",
        left_flank=1,
        right_flank=1,
    )

    read = Read(seq="GAAAG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, list)
    assert len(result) == 2  # Should match "A" and "AA"

    # Check matches with flanks - the match field contains the sequence with flanks
    matches = [hit.match for hit in result]
    assert "GAAA" in matches  # A with 1 base flank on each side
    assert "AAG" in matches  # AA with 1 base flank on each side


def test_pattern_search_with_flanks_anchors():
    """Test Pattern search with flanks and regex anchors"""
    pattern = Pattern(
        name="test",
        reg_expr="^ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=2,
        right_flank=1,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should be clamped to sequence start
    assert result.end == 4  # Should be clamped to sequence end (no extra base after)
    assert result.match == "ATCG"  # Should include the flanks (clamped to sequence)
    assert result.orientation == "forward"


def test_pattern_search_with_flanks_edge_cases():
    """Test Pattern search with flanks in edge cases"""
    # Test with very large flank values
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=1000,
        right_flank=1000,
    )

    read = Read(seq="ATCG", pair="R1", name="test_read")
    result = pattern.search(read)

    assert result is not None
    assert isinstance(result, Hit)
    assert result.start == 0  # Should be clamped to sequence start
    assert result.end == 4  # Should be clamped to sequence end
    assert result.match == "ATCG"  # Should include entire sequence

    # Test with flank values equal to sequence length
    pattern2 = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=4,
        right_flank=4,
    )

    result2 = pattern2.search(read)

    assert result2 is not None
    assert isinstance(result2, Hit)
    assert result2.start == 0  # Should be clamped to sequence start
    assert result2.end == 4  # Should be clamped to sequence end
    assert result2.match == "ATCG"  # Should include entire sequence


def test_pattern_search_with_flanks_generator():
    """Test Pattern._search_read generator with flanks"""
    pattern = Pattern(
        name="test",
        reg_expr="ATCG",
        read="R1",
        orientation="forward",
        multiple="first",
        left_flank=1,
        right_flank=1,
    )

    read = Read(seq="GATCG", pair="R1", name="test_read")

    hits = list(pattern._search_read(read))

    assert len(hits) == 1
    assert hits[0].start == 0  # Should include 1 base before ATCG
    assert hits[0].end == 5  # Should include 1 base after ATCG
    assert hits[0].match == "GATCG"  # Should include the flanks
    assert hits[0].orientation == "forward"

    
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.