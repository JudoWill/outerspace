"""configuration for defining motifs"""
from tomlkit import comment
from tomlkit import document
from tomlkit import nl
from tomlkit import table
from tomlkit import toml_file

print(f'CCCCCCCCC config({__name__})')

# Format taken from TOML Quickstart Kit
class Cfg:
    """configuration for defining motifs"""
    
    # Upon Cfg() __init__() runs - bc we are creating document here
    def __init__(self):
        self.doc = self._init_doc()

    def __str__(self):
        return self.doc.as_string()

    @staticmethod
    def _init_doc():
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
        
        grna.add("regex_flags", "BESTMATCH")
        grna.add("umi_pattern_forward", '(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}') 
        grna.add("protospacer_pattern", '(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')
        
        grna.add("back_umi_forward", 'gtgtgtcagttagggtgtggaa')
        
        grna.add("umi_pattern_reverse",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        
        # Adding the table to the document
        doc.add("define_motifs", grna)
        return doc
