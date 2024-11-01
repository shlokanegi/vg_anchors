class GafReader:
    """
    This class processes the gaf alignment file from Giraffe HiFi and return the necessary data
    structured for "fast" anchor queries.
    It takes the start/end of the alignment in the path, the path and the cs:Z optional tag.
    It parses the path field to record nodes and orientations, stored in two numpy arrays (int and bool)
    For gaf tags see: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
    For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
    """
    
    def __init__(self, text_file_path: str):
        self.file_path = text_file_path

    def get_lines(self):
        with open(self.file_path, 'r') as f:
            for line in f:
                yield line.strip()
