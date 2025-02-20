"""configuration for defining motifs"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "??"

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
    def __init__(self, filename=None):
        self.filename = filename

    #file does not need to be given name
    def read_file(self):
        """Read the file specified"""
        print(f'  READ: {self.filename}')
        return TOMLFile(self.filename).read() if filename is not None else None


    def write_file(self):
        """Write a default configuration file"""
        doc = self.get_doc_default()
        if self.filename is not None:
            TOMLFile(self.filename).write(doc)
            print(f'  WROTE: {self.filename}')
        else:
            print(f'  PLEASE PROVIDE FILE NAME IN ORDER TO WRITE')
        return doc

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
### Next need to finish


    def get_doc_default(self):
        """Getting default config document object"""
        doc = document()
        doc.add(comment("Configuration file for gRNA extraction."))
        doc.add(nl())
        doc.add("title", "gRNA Extraction Input")
        doc.add(nl())
        # Using doc["title"] = "TOML Example" is also possible

        owner = table()
        owner.add("name", "SCB, DK, WND, RB")
        owner.add("organization", "DUCOM")
        # Adding the table to the document
        doc.add("owner", owner)
        doc.add(nl())

        self._addregxlist(doc)
        doc.add(nl())
        self._addregxlist1(doc)
        doc.add(nl())
        self._addregxlist2(doc)
        doc.add(nl())

        # This is made into an array
        #regxlist = table()
        doc.add(comment("Edit this list to specify what part of sequence you want to capture."))
        doc.add(comment("Named patterns/groups captured  will be saved in a .csv file. Ex: (?P<UMI>.{8})"))

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

    #regxlist here applies to pattern in both reads
    @staticmethod
    def _addregxlist(doc):
        """regxlist here applies to pattern in both reads"""
        arr = array()
        arr.multiline(True)
        doc['regxlist'] = arr

    #regxlist here applies to pattern in first read of paired read
    @staticmethod
    def _addregxlist1(doc):
        """regxlist1 here applies to pattern in first read of paired read"""
        arr = array()
        arr.multiline(True)
        doc['regxlist1'] = arr

        #regxlist.add("regex_flags", "BESTMATCH")

        doc['regxlist1'].add_line('(?P<UMI>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}')

        doc['regxlist1'].add_line('(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}'
                                 '(?P<protospacer>.{19,21})'
                                 '(?P<protospacer2>GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}')

    @staticmethod
    def _addregxlist2(doc):
        """regxlist2 here applies to pattern in 2nd read of paired read"""
        arr = array()
        arr.multiline(True)
        doc['regxlist2'] = arr

    def __str__(self):
        return self.doc.as_string()



 # Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
