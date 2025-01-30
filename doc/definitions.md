# Definitions
```
forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}'    , flags=regex.BESTMATCH)
protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protosp    acer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}', flags=regex.BESTMATCH)
back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
back_umi_rc = reverse_complement(back_umi_forward)
reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', fl    ags=regex.BESTMATCH)
```

