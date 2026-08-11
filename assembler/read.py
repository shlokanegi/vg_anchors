from sys import stderr
from assembler.anchor import Anchor

class Read:

    def __init__(self, name: str, strand: str) -> None:
        self.name = name
        self.strand = strand
        self.journey = []   # anchor objects in the order they were visited by the read

    def add_anchor(self, anchor: Anchor) -> None:
        self.journey.append(anchor)
        
    
