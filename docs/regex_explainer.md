# Regex Explanation
  
A `regex` or regular expression is a sequence of characters that defines a *pattern* that can be searched across input sequences.
Regular expressions can be complex, but only a small fraction of the language is needed for our purposes.

 - Constants : Any letter (A, C, G, T) is matched as is.
 - `[AGC]` : Square brackets indicate that any contained characters will match. Useful for ambiguity.
 - `.` : A dot character matches any letter.
 - `(AGTA)` : Group characters so subsequent modifiers apply to multiple characters
 - `{n}`, `{a,b}` : Matches the previous character/group exactly `n` times or at least `a` times but not more than `b` times.
 - `?P<name>` : Names a group for later extraction.
 - `{s<=a}` : Defines fuzzy matching of the previous group. `s` - substitutions, `i` - insertions, `d` - deletions, `e` - errors.

These will account for the majority of pieces needed to create your pattern.

## Problem Definition

Imagine you are performing an NGS experiment in which you used PCR primers tagged with unique molecular indices (UMIs).
This will tag each first-round PCR product with a unique index on each end of the molecule.

Considering only the forward primer:

### Capture

Forward:
```
agtacgtacgtagctagNNSNNSNNSagtacgtacgtacgatttagctagtacg
```

There are many ways one could design a regex which extracts the UMI from this technical sequence.
The simplest is to use the constant regions surrounding the tag.

Like so:
```
gtagctag.{9}agtacgta
```

or 
```
gtagctag(..[ACG]){3}agtacgta
```

### Named groups

The pattern pattern above will match the the constant sequence on either side as well as the UMI.
You can define which part of your pattern is relevant using named groups.

To do this, we encase our target as a 'group' and then add `?P<name>` to the front.

```
gtagctag(?P<UMI_F>(..[ACG]){3})agtacgta
```
By default, `outerspace` will extract all named groups as individual entities.
This is useful if you have multiple patterns to match in a single read.

### Fuzzy Matching

There are certain adapter designs and sequencing modalities where it is useful to allow for an approximate match.
With fuzzy matching, you can create patterns which tolerate sequencing or biological variability. 
These can be defined at the group-level allowing you to set different tolerances at different parts of the pattern.

Imagine the experiment was performed on a nanopore sequencer with with the possibility of mismatches and indels introduced during sequencing.
This could be accounted for with fuzzy matching:

```
(gtagctag){e<=2}(?P<UMI_F>(..[ACG]){3})(agtacgta){e<=2}
```

This will allow for up-to two errors (mismatch, insertion, or deletion) during the search.
Beware, searching with mismatches drastically increases the complexity and slows down the search.


Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.