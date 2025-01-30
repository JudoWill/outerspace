# gRNA EXTRACTION PROJECT  
  
## Description:
This pipeline is intended to to extract the gRNA protospacer from reads. In future will adjust to also extract tracrRNA from reads for a CRISPR Screen.  
Original Location: Mistake-Not - /share/nonn-lab/Rachel-test-cripr 

## Tools utilized
### Regex  
To find the protospacer because it contains a surrounding pattern.  
- Regex Link: -	https://pypi.org/project/regex/ 
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
    - (?:foo){i} match “foo”, permitting insertions
    - (?:foo){d} match “foo”, permitting deletions
    - (?:foo){s} match “foo”, permitting substitutions
    - (?:foo){i,s} match “foo”, permitting insertions and substitutions
    - (?:foo){e} match “foo”, permitting errors
- Additional parameters can be added - check on website
    
###     



