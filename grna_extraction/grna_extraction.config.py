from tomlkit import comment
from tomlkit import document
from tomlkit import nl
from tomlkit import table

# Format taken from TOML Quickstart Kit

doc = document()
doc.add(comment("TOML document for gRNA extraction."))
doc.add(nl())
doc.add("title", "gRNA Extraction Input")
# Using doc["title"] = "TOML Example" is also possible

owner = table()
owner.add("name", "SCB, DK, WND, RB")
owner.add("organization", "DUCOM")

# Adding the table to the document
doc.add("owner", owner)

grna = table()
grna.add("name", "SCB, DK, WND, RB")
grna.add("organization", "DUCOM")

grna.add("forward_umi_reg = regex.compile('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}', flags=regex.BESTMATCH) 
grna.add("protospacer_reg = regex.compile('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}', flags=regex.BESTMATCH)
grna.add("back_umi_forward = 'gtgtgtcagttagggtgtggaa'.upper()
grna.add("back_umi_rc = reverse_complement(back_umi_forward)
grna.add("reverse_umi_reg = regex.compile(f'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}', flags=regex.BESTMATCH) 
# Adding the table to the document
doc.add("define_motifs", grna)
