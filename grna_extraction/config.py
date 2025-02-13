"""configuration for defining motifs"""
from tomlkit import comment
from tomlkit import document
from tomlkit import nl
from tomlkit import table
from tomlkit import array
from tomlkit.toml_file import TOMLFile



print(f'CCCCCCCCC config({__name__})')

# Format taken from TOML Quickstart Kit
class Cfg:
    """configuration for defining motifs"""
    
    # Upon Cfg() __init__() runs - bc we are creating document here
    # referring to self and filename in arguments that are passed by the user
    # Need self in the class to have acces to everything within the class
    def __init__(self, filename):
        self.filename = filename

    def read_file(self):
        """Read the file specified"""
        print(f'  READ: {self.filename}')
        return TOMLFile(self.filename).read()
        

    def write_file(self):
        """Write a default configuration file"""
        doc = self.get_doc_default()
        TOMLFile(self.filename).write(doc)
        print(f'  WROTE: {self.filename}')

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
        """Get forward protospacer pattern"""
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
        doc.add(nl())
        # Using doc["title"] = "TOML Example" is also possible
        
        owner = table()
        owner.add("name", "SCB, DK, WND, RB")
        owner.add("organization", "DUCOM")
        # Adding the table to the document
        doc.add("owner", owner)
        
        #regxlist = table()
        doc.add(comment("Edit this list to specify what part of sequence you want to capture."))
        doc.add(comment("Named patterns/groups captured  will be saved in a .csv file. Ex: (?P<UMI>.{8})"))
        arr = array()
        arr.multiline(True)
        doc['regxlist'] = arr
        
        #regxlist.add("regex_flags", "BESTMATCH")
        
        #regxlist.add("umi_pattern_forward", '(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}') 
        doc['regxlist'].add_line('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}') 
        # pulled apart this variable into 3 parts
        #regxlist.add("umi_pattern_forward_num_umi_nt", 8) 
        #regxlist.add("umi_pattern_forward_pattern_nt", 'CTTGGCTTTATATATCTTGTGG') 
        #regxlist.add("umi_pattern_forward_mismatch_max_nt", 4) 


        #regxlist.add("protospacer_pattern", '(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')
        doc['regxlist'].add_line('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?:GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')
        # pulled apart this variable into 6 parts
        # NOTE THAT THE PROTOSPACER ITSLEF CAN BE A RANGE OR JUST ONE # - WILL NEED TO MAKE THIS SO THAT IT TESTS APPROPRIATELY
        #regxlist.add("pattern_forward_upstream_protospacer", 'TATCTTGTGGAAAGGACGAAACACC')
        #regxlist.add("pattern_forward_upstream_protospacer_mismatch_max_nt", 4)
        #regxlist.add("pattern_forward_num_protospacer_nt_range_from", 19)
        #regxlist.add("pattern_forward_num_protospacer_nt_range_to", 21)
        #regxlist.add("pattern_forward_downstream_protospacer", 'GTTTAAGTACTCTGTGCTGGAAACAG')
        #regxlist.add("pattern_forward_protospacer_mismatch_max_nt", 4)
        

         
        #regxlist.add("back_umi_forward", 'gtgtgtcagttagggtgtggaa')
        #regxlist.add('gtgtgtcagttagggtgtggaa')
        # Kept this variable in 1 part
        # regxlist.add("umi_pattern_forward_downstream_nt", 'gtgtgtcagttagggtgtggaa')

##############FINISH THIS LAST ONE BELOW & MAKE SURE TO MAKE CHANGES AND ADD TESTS TO TESTS SCRIPT & MAKE SURE TO ADD TO THE DICTIONARY IN THAT FILE AS WELL#############
        
        #regxlist.add("umi_pattern_reverse",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        #regxlist.add('(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        # Kept this variable in  3 parts
        #regxlist.add("umi_pattern_reverse_downstream_nt",8) 
        #regxlist.add("umi_pattern_reverse",'(?P<UMI>.{{8}})(?:{back_umi_rc}){{s<=4}}') 
        
        # Adding the table to the document
        #doc.add("define_motifs", regxlist)
        return doc

    def __str__(self):
        return self.doc.as_string()
