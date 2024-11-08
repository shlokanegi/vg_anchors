import numpy as np


class GafReader:
    """
    This class processes the gaf alignment file from Giraffe HiFi and return the necessary data
    structured for "fast" anchor queries.
    It takes the start/end of the alignment in the path, the path and the cs:Z optional tag.
    It parses the path field to record nodes and orientations, stored in two numpy arrays (int and bool)
    It parses the cs tag to represent it as a list of steps. Each step is a tuple representing the
    difference operator (+,-,:,=,*) and the length of the operation in basepairs
    """

    def __init__(
        self,
        read_id: str,
        nodes_array: np.ndarray,
        orientation_array: np.ndarray,
        path_start: int,
        path_end: int,
        relative_strand: bool,
        cs_string: str,
    ):
        """
        Constructor for initializing the path object with the needed tags
        """
        self.read_id: int = read_id
        self.nodes_array = nodes_array
        self.orientation_array = orientation_array
        self.path_start = path_start
        self.path_end = path_end
        self.relative_strand = relative_strand
        self.cs_string = cs_string

    def __repr__(self):
        """
        A string representation of the object showing some info about the path.
        """
        strand: str = "+" if self.relative_strand else "-"
        return (
            f"Read name: {self.read_id}, "
            f"Path start (node_id): {self.path_start})"
            f"Path end (node_id): {self.path_end})"
            f"Relative strand: {strand}"
        )


def gaf_sanity_check(gaf_line: str) -> bool:
    """
    This function returns True if the gaf line is consistent to the one expected by the class.
    It cannot verify that the nodes are in the gfa format.
    """
    line_elements = gaf_line.strip().split("\t")
    if len(line_elements) != 16:
        return False

    # there should be the 'cs' optional field tag
    if line_elements[15].startswith("cs:Z:"):
        return True

    return False


def get_path_from_line(gaf_line: str):
    """
    This function returns a properly initialized Path object
    """
    line_elements = gaf_line.strip().split()
    read_name = line_elements[0]
    relative_strand = True if line_elements[4] == "+" else False
    path_start = line_elements[7]
    path_end = line_elements[8]
    map_q = line_elements[11]

    nodes_list = []
    orientation_list = []
    for i in range(0, len(line_elements[5]), 2):
        orientation_list.append(line_elements[5][i])
        nodes_list.append(line_elements[5][i + 1])

    cs_line = [i for i in parse_cs_line(line_elements[15])]

    for el in cs_line:
        print(f"Flag: {el[0]} ; length {el[1]}")

    # read_name, nodes_list, orientations_list, path_start, path_end, strandness
    # return (read_name, np.array(nodes_list, dtype=np.int32), np.array(orientations_list, dtype=np.bool_), path_start, path_end, strandness)
    return Path()


def parse_cs_line(cs_string: str):
    flag_chars = ":*+-="
    i = 5

    while i != len(cs_string):
        if cs_string[i] in "+-=":
            flag_c = cs_string[i]
            count = i
            i += 1
            while i < len(cs_string) and cs_string[i] not in flag_chars:
                i += 1
            yield (flag_c, i - count - 1)
            continue

        elif cs_string[i] == ":":
            number = ""
            i += 1
            while cs_string[i] not in flag_chars:
                if i < len(cs_string):
                    number += cs_string[i]
                    i += 1
            yield (":", int(number))
            continue

        elif cs_string[i] == "*":
            yield ("*", 1)
            i = i + 3
            continue
