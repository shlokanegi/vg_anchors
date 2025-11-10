
class Node:

    def __init__(self,id,length,orientation) -> None:
        self.id: int = id
        self.length: int = length
        self.orientation: bool = orientation
        # When we extend inside try_extension, when we try to extend right, if node orientation is False, 
        # we have to actually pass extend_left = True to get_degree() and follow_edges()

    def __eq__(self, other_node) -> bool:
        """
        Checks equality of node objects based on Node IDs only
        """

        return (True if self.id == other_node.id else False)
