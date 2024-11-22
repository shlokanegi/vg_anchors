from sys import stderr
class Anchor:

    def __init__(self) -> None:
        self._nodes: list = []
        self.snarl_id: int = 0
        self.genomic_position: int = 0
        self.baseparilength: int = 0
        self.num_sequences: int = 0

    def __len__(self):
        return len(self._nodes)
    
    def __getitem__(self, position):
        return self._nodes[position]

    def add(self, node):
        self._nodes.append(node)

    def add_snarl_id(self, snarl_id) -> None:
        self.snarl_id = snarl_id
    
    def __repr__(self) -> str:
        anchor_str = ""
        for node in self._nodes:
            orientaiton = ">" if node.orientation else "<"
            anchor_str += orientaiton + str(node.id)
        return anchor_str

    def bandage_representation(self) -> str:
        bandage_nodes_str = ""
        for node in self._nodes:
            bandage_nodes_str += "," + str(node.id)
        return bandage_nodes_str[1:]
    
    def add_sequence(self) -> None:
        self.num_sequences += 1

    def compute_bp_length(self) -> int:
        """
        This function computes the size, in basepairs, of an anchor.

        Returns
        -------
        anchor_size: int
        The length in basepairs of the anchor
        """
        anchor_str = ""
        for node in self._nodes:
            orientaiton = ">" if node.orientation else "<"
            anchor_str += orientaiton + str(node.id)
        self.baseparilength = (self._nodes[0].length // 2) + (self._nodes[-1].length // 2)
        print(f"Anchor {anchor_str}: Start_bp: {(self._nodes[0].length // 2)}, end_bp: {(self._nodes[-1].length // 2)}",file=stderr, end= " ")

        for node_handle in self._nodes[1:-1]:
            self.baseparilength += node_handle.length
            print(f"+ {node_handle.length}",end=" ", file=stderr)
        print(f"total: {self.baseparilength}",file=stderr)
        return 

    def get_sentinel_id(self) -> int:
        """
        This function takes a path traversing a snarl (candidate anchor) and returns the sentinel associated to it. The sentinel is not dependent on the anchor orientation.

        Returns
        -------
        sentinel: int
        the node_id of the sentinel node
        """
        # Determine the index based on the orientation of the first element
        if self._nodes[0].orientation:
            # If reversed, count from the end
            index = -((len(self._nodes) + 1) // 2)
        else:
            # If not reversed, count from the beginning
            index = (len(self._nodes) - 1) // 2

        return self._nodes[index].id

    def __eq__(self, other_anchor) -> bool:
        """
        This function verifies if two paths in the same snarls are equal, independently of the orientation.

        Parameters
        ----------
        path_1: list
        list of node handle objects of the first path
        path_2: list
        list of node handle objects of the first path

        Returns
        -------
        bool
        True if the paths are equal else False
        """

        # if different in length, don't need to bother
        if len(self._nodes) != len(other_anchor):
            return False

        pos_1 = pos_2 = 0
        orientation_concordance = self._nodes[pos_1].orientation == other_anchor[pos_2].orientation
        if not orientation_concordance:
            pos_1 = len(self._nodes) - 1

        for _ in range(len(self._nodes)):
            if self._nodes[pos_1].id != other_anchor[pos_2].id or orientation_concordance != (
                self._nodes[pos_1].orientation == other_anchor[pos_2].orientation
            ):
                return False
            pos_1 += 1 if orientation_concordance else -1
            pos_2 += 1
        return True