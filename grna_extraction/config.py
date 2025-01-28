"""configuration for defining motifs"""
from tomlkit import comment
from tomlkit import document
from tomlkit import nl
from tomlkit import table

print(f'CCCCCCCCC config({__name__})')

# Format taken from TOML Quickstart Kit
class Cfg:
    """configuration for defining motifs"""

    @staticmethod
    def get_doc():
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
        
        grna.add("regex_flags", "BESTMATCH")
        grna.add("forward_umi_pattern", '(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}') 
        grna.add("protospacer_pattern", '(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')
        
        grna.add("back_umi_forward", 'gtgtgtcagttagggtgtggaa')
        
        grna.add("reverse_umi_pattern",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        
        # Adding the table to the document
        doc.add("define_motifs", grna)
