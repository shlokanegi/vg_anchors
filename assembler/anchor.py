from sys import stderr
class Anchor:

    def __init__(self) -> None:
        self._nodes: list = []
        self.snarl_id: int = 0
        self.genomic_position: int = 0
        self.baseparilength: int = 0
        self.num_sequences: int = 0
        self.chromosome: str = ""
        self.reference_paths_covered: list = []
        # self.path_matched_reads: list = []
        self.bp_matched_reads: list = []

    def add_reference_path(self, path):
        self.reference_paths_covered.append(path)

    def __len__(self):
        return len(self._nodes)
    
    def __getitem__(self, position):
        return self._nodes[position]

    def add(self, node):
        self._nodes.append(node)

    def merge_anchor(self, new_anchor, insert_left=False) -> bool:
        if insert_left:
            if self._nodes[0].id != new_anchor[-1].id:
                return False
            self._nodes[:0] = new_anchor[:-1]
        else:
            if self._nodes[-1].id != new_anchor[0].id:
                return False
            self._nodes.extend(new_anchor[1:])
        return True

    def flip_anchor(self):
        self._nodes = self._nodes[::-1]
        for node in self._nodes:
            node.orientation = not node.orientation

    # def insert(self, node, index):
    #     self._nodes.insert(index, node)

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
        First and last node (the ones at the boundaries) are just half the node, while the
        nodes in the middle are taken in full.

        Returns
        -------
        anchor_size: int
        The length in basepairs of the anchor
        """

        self.baseparilength = (self._nodes[0].length // 2) + (self._nodes[-1].length // 2)
        for node_handle in self._nodes[1:-1]:
            self.baseparilength += node_handle.length
        # self.baseparilength = sum([node_handle.length for node_handle in self._nodes])
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
    
    def get_reference_paths(self):
        out_s = ""
        for el in self.reference_paths_covered:
            out_s += el + ','
        return out_s[:-1]  

    def get_bed(self):
        #CHROM CHROM_START CHROM_END NAME
        return f"{self.chromosome}\t{self.genomic_position}\t{self.genomic_position+self.baseparilength}\t{self.__repr__}"