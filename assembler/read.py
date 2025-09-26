from sys import stderr
from assembler.anchor import Anchor

class Read:

    def __init__(self, name: str) -> None:
        self.name = name
        self.snarls_to_anchors = {}   # {snarl_id: [anchor1, anchor2, ...]}
    
    def add_anchor(self, anchor):
        """
        Adds anchor to the snarl_to_anchors dictionary of the read object.
        """
        snarl_id = anchor.snarl_id

        if snarl_id not in self.snarls_to_anchors:
            self.snarls_to_anchors[snarl_id] = []
        self.snarls_to_anchors[snarl_id].append(anchor)
