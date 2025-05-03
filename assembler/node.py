
class Node:

    def __init__(self,id,length,orientation) -> None:
        self.id: int = id
        self.length: int = length
        self.orientation: bool = orientation

    def __eq__(self, other_node) -> bool:
        """
        Checks equality of node objects based on Node IDs only
        """

        return (True if self.id == other_node.id else False)
