# Regex Explanation
  
A regex or regular expression is a sequence of characters that defines a pattern that can be utilized in this tool to find and capture exact match when looking for a precise pattern(s), or "fuzzy" match which allows for error(s). 

### start with regular stuff

### Example  
```
Regex: (?P<planet>.{8})(?:CTTGGCTAA){s<=3}  
Read Sequence: N N N N N N N N N N N N N N N N N N C T C G G C G T A N N N N N  
                                  |---------------|-----------------|  
                                   (?P<planet>.{8}) (?:CTTGGC){s<=3}
  
If match is successful: planet = NNNNNNNN  
```
### Capture name
```
(?P<name>...)
```
Description: Anything identified within the parenthesis is matched and captured with the name you provide.
  
### Capture exact # of characters and name them
```
(?P<name>.{8})
```
Description: Capture any 8 characters and name group.
  
### Non-capturing group
```
(?:...)
```
Description: If don't need to capture and just want to match.
```
(?:CTTGGC)
```
Example: Match exact CTTGGC sequence with no capture.

### Match with errors allowed
```
(?:CTTGGC){s<=3}
```
Description: Match CTTGGCT sequence allowing up to 3 substitutions.
  
### Match sequence with errors allowed and capture exact # of characters to a named group
```
(?P<name>.{8})(?:CTTGGC){s<=3}
```
Description: Capture first 8 characters and name group. Then match TTGGC and allow up to 3 substitutions. 
  
If match is successful name = NNNNNNNN


  
### Walkthrough Example
Example described from [walkthrough](http://outerspace/docs/walkthrough.md)     
```
Regex: (?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}    
Description: Capture first 8 characters and name UMI_5prime. Then match CTTGGCTTTATATATCTTGTGG and allow up to 4 substitutions.
``` 
  
```                                                             
Regex: (?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})  
Description: Match TATCTTGTGGAAAGGACGAAACACC allowing up to 4 substitutions. Then capture following characters with a length of 19-21 and name that group protospacer.     
```
  
```    
Regex:(?P<UMI_3prime>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}  
Description: Capture first 8 characters and name group UMI_3prime. Then match TTCCACACCCTAACTGACACAC allowing up to 4 substitutions.
```

### Regex
- [Regex Link](https://pypi.org/project/regex/) 
- Regex has approximate fuzzy matching
- Regex usually attempts an exact match, but sometimes an approximate, or “fuzzy”, match is needed, for those cases where the text being searched may contain errors in the form of inserted, deleted or substituted characters.

- A fuzzy regex specifies which types of errors are permitted, and, optionally, either the minimum and maximum or only the maximum permitted number of each type. (You cannot specify only a minimum.)

- The 3 types of error are:
    - Insertion, indicated by “i”
    - Deletion, indicated by “d”
    - Substitution, indicated by “s”
- In addition, “e” indicates any type of error.
- The fuzziness of a regex item is specified between “{” and “}” after the item.
- Examples:
    - foo match “foo” exactly
    - (?:AAA){i} match “AAA”, permitting insertions
    - (?:AAA){d} match “AAA”, permitting deletions
    - (?:AAA){s} match “AAA”, permitting substitutions
    - (?:AAA){i,s} match “AAA”, permitting insertions and substitutions
    - (?:AAA){e} match “AAA”, permitting errors
- Additional parameters can be added - check on website



### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.