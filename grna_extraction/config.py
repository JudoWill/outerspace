"""configuration for defining motifs"""
from tomlkit import comment
from tomlkit import document
from tomlkit import nl
from tomlkit import table
from tomlkit.toml_file import TOMLFile


print(f'CCCCCCCCC config({__name__})')

# Format taken from TOML Quickstart Kit
class Cfg:
    """configuration for defining motifs"""
    
    # Upon Cfg() __init__() runs - bc we are creating document here
    #def __init__(self):
    #    self.doc = self.get_doc_default()

    def write_file(self, filename):
        """Write a default configuration file"""
        doc = self.get_doc_default()
        TOMLFile(filename).write(doc)
        print(f'  WROTE: {filename}')

    def get_umi_pattern_forward(self):
        """Get forward UMI regex pattern"""
        return ('(?P<UMI>.{'
            f"{self.doc['define_motifs']['umi_pattern_forward_num_umi_nt']}"
            '})(?:'
            f"{self.doc['define_motifs']['umi_pattern_forward_pattern_nt']}"
            '){s<='
            f"{self.doc['define_motifs']['umi_pattern_forward_mismatch_max_nt']}"
            '}')

    def get_protospacer_forward(self):
        """Get forward prtospacer pattern"""
        return ('(?:'
            f"{self.doc['define_motifs']['pattern_forward_upstream_protospacer']}"
            '){s<='
            f"{self.doc['define_motifs']['pattern_forward_upstream_protospacer_mismatch_max_nt']}"
            '}(?P<protospacer>.{'
            f"{self.doc['define_motifs']['pattern_forward_num_protospacer_nt_range_from']}"
            ','
            f"{self.doc['define_motifs']['pattern_forward_num_protospacer_nt_range_to']}"
            '})(?:'
            f"{self.doc['define_motifs']['pattern_forward_downstream_protospacer']}"
            '){s<='
            f"{self.doc['define_motifs']['pattern_forward_protospacer_mismatch_max_nt']}"
            '}')
        
    def get_umi_pattern_forward_downstream_nt(self):
        """Get downstream forward UMI regex pattern"""
        return f"{self.doc['define_motifs']['umi_pattern_forward_downstream_nt']}"


##########STOPPED HERE#####
### Need to fix this error bc when run tests script its giving string error
### Next need to finish 


    @staticmethod
    def get_doc_default():
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
        
        #grna.add("regex_flags", "BESTMATCH")
        
        grna.add("umi_pattern_forward", '(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}') 
        # pulled apart this variable into 3 parts
        #grna.add("umi_pattern_forward_num_umi_nt", 8) 
        #grna.add("umi_pattern_forward_pattern_nt", 'CTTGGCTTTATATATCTTGTGG') 
        #grna.add("umi_pattern_forward_mismatch_max_nt", 4) 


        grna.add("protospacer_pattern", '(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')
        # pulled apart this variable into 6 parts
        # NOTE THAT THE PROTOSPACER ITSLEF CAN BE A RANGE OR JUST ONE # - WILL NEED TO MAKE THIS SO THAT IT TESTS APPROPRIATELY
        #grna.add("pattern_forward_upstream_protospacer", 'TATCTTGTGGAAAGGACGAAACACC')
        #grna.add("pattern_forward_upstream_protospacer_mismatch_max_nt", 4)
        #grna.add("pattern_forward_num_protospacer_nt_range_from", 19)
        #grna.add("pattern_forward_num_protospacer_nt_range_to", 21)
        #grna.add("pattern_forward_downstream_protospacer", 'GTTTAAGTACTCTGTGCTGGAAACAG')
        #grna.add("pattern_forward_protospacer_mismatch_max_nt", 4)
        

         
        grna.add("back_umi_forward", 'gtgtgtcagttagggtgtggaa')
        # Kept this variable in 1 part
        # grna.add("umi_pattern_forward_downstream_nt", 'gtgtgtcagttagggtgtggaa')

##############FINISH THIS LAST ONE BELOW & MAKE SURE TO MAKE CHANGES AND ADD TESTS TO TESTS SCRIPT & MAKE SURE TO ADD TO THE DICTIONARY IN THAT FILE AS WELL#############
        
        grna.add("umi_pattern_reverse",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        # Kept this variable in  3 parts
        #grna.add("umi_pattern_reverse_downstream_nt",8) 
        #grna.add("umi_pattern_reverse",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        
        # Adding the table to the document
        doc.add("define_motifs", grna)
        return doc

    def __str__(self):
        return self.doc.as_string()
