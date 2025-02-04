# Definitions
### NEED TO ADD MORE DETAILS ABOUT THE CODE ITSELF
```
forward_umi_reg
    = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}    , flags=regex.BESTMATCH)`
    = Primer (Read 1 SEQ Primer) complementary to U6 promoter for gRNA overlapped by forward UMI (LIB gRNA query UMI FWD) allowing 4 substitutions 

protospacer_reg 
    = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protosp    acer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}', flags=regex.BESTMATCH)
    = Obtaining protospacer sequence
    = Protospacer the length of 19-21
    = Sequences surrounding protospacer: 
        = Upstream Sequence within U6 promoter for gRNA, partial overlap of primer (Read 1 SEQ Primer) and partial overlap of forward UMI (LIB gRNA query UMI FWD)
        = Downstream Sequence within the SaCas9 gRNA tracrRNA/ scaffold

back_umi_forward
    = 'gtgtgtcagttagggtgtggaa'.upper()
    = Sequence within the start of SV40 promoter for Cas9 (downstream the gRNA), complementary to the primer (Read2 SEQ Primer), overlapping reverse UMI (LIB gRNA query UMI REV)

back_umi_rc = reverse_complement(back_umi_forward)
    = reverse complement the sequence within the start of SV40 promoter

reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', fl    ags=regex.BESTMATCH)
    = allow 4 substitutions for the reverse complement of sequence withing start of SV40 promoter
```

