# data_processor.py
import numpy as np

class Path:
    """
    This class encapsulates two arrays: one for integers and one for signs (+/-),
    and provides methods for internal manipulations of the data.
    """
    
    def __init__(self, read_id: str, nodes_array: np.ndarray, orientation_array: np.ndarray, path_start: int, path_end: int, relative_strand: bool):
        """
        Constructor for initializing the path object with the needed info.
        """
        self.read_id = read_id
        self.nodes_array = nodes_array
        self.orientation_array = orientation_array
        self.path_start = path_start
        self.path_end = path_end
        self.relative_strand = relative_strand

    def get_range(self) -> np.ndarray:
        """
        Returns the range of the alignment in the graph as a tuple (path_start,path_end)
        """
        return (self.path_start,self.path_end)
    
    def __repr__(self):
        """
        A string representation of the object showing some info about the path.
        """
        strand:str = "+" if self.relative_strand else "-"
        return (f"Read name: {self.read_id}, "
                f"Path start (node_id): {self.path_start})"
                f"Path end (node_id): {self.path_end})"
                f"Relative strand: {strand}"
                )

# If you want to define helper functions here, you can too.
# For example:
def get_path_from_line(read_name, nodes_list, orientations_list, path_start, path_end, strandness):
    """
    Helper function to convert basic lists into NumPy arrays and
    return a DataProcessor example.
    """
    return Path(read_name, np.array(nodes_list, dtype=np.int32), np.array(orientations_list, dtype=np.bool_), path_start, path_end, strandness)