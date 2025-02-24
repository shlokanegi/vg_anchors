from sys import stderr
from assembler.constants import (
    EXPECTED_GAF_TAGS,
    EXPECTED_MAP_Q,
    MIN_CS_LEN,
    MAP_Q_ID,
    READ_NAME_ID,
    RELATIVE_STRAND_ID,
    READ_LEN,
    PATH_START_ID,
    PATH_END_ID,
    PATH_ID,
    CS_TAG_ID,
)

"""
This functions process the gaf alignment file from Giraffe HiFi and return the necessary data
structured for "fast" anchor queries.
- process_line takes start/end of the alignment in the path, the path and the cs:Z optional tag.
It parses the path field to record nodes and orientations, stored in two numpy arrays (int and bool)

- parse_cs_line parses the cs:Z field to structure the cigar information to verify the basepair alignment between the path sequence and the read 

For gaf tags see: https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf
For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
"""


def processGafLine(gaf_line: str):
    """
    It parses a GAF line to extract and structure useful tags and returns them in a list

    Parameters
    ----------
    gaf_line : string
        a gaf line stripped of whitespaces at the beginning and end

    Returns
    -------
    list
        list of processed tags:
        read_name : string
        read_len : int
        relative_strand : bool (True if +, False else)
        path_start : int - start of the alignment in the path sequence
        path_end : int - end of the alignment in the path sequence
        nodes_list : list - of node_ids of the nodes walked by the path
        orientation_list : list - of node orientations of the nodes walked by the path
        cs_line : list - succession of tuples describing the cigar
    """

    line_elements = gaf_line.split()
    #print(line_elements[1:2])
    if not(line_elements[1].isnumeric()):
        #print('Not numeric')
        del line_elements[1:3]
    #print(f"# el: {len(line_elements)}, expected_tags: {EXPECTED_GAF_TAGS}")
    # First verify that the gaf line contains an usable alignment
    if (len(line_elements) == EXPECTED_GAF_TAGS) and int(
        line_elements[MAP_Q_ID]
    ) >= EXPECTED_MAP_Q:
        #print("processing")
        # extract needed tags
        read_name = line_elements[READ_NAME_ID]
        read_len = int(line_elements[READ_LEN])
        #relative_strand = True if line_elements[RELATIVE_STRAND_ID] == "+" else False
        path_start = int(line_elements[PATH_START_ID])
        path_end = int(line_elements[PATH_END_ID])
        
        # decompose the path into nodes and orientations arrays
        nodes_list = []
        orientation_list = []

        curr_node_string = ""
        for char in line_elements[PATH_ID]:
            if char in "><":
                orientation_list.append(True if char == ">" else False)
                if len(curr_node_string) != 0:
                    nodes_list.append(int(curr_node_string))
                    curr_node_string = ""
            else:
                curr_node_string += char
        if len(curr_node_string) != 0:
            nodes_list.append(int(curr_node_string))
        
        count_positive_orientation_nodes = orientation_list.count(True)
        relative_strand = True if count_positive_orientation_nodes > (len(orientation_list) / 2) else False

        # decompose the cs tag into alignment steps
        if len(line_elements[CS_TAG_ID]) > MIN_CS_LEN:
            cs_line = [i for i in parse_cs_tag(line_elements[CS_TAG_ID])]
        else:
            print("ERROR IN CS LINE.",flush=True, file=stderr)
            return None
        # print(f"{read_name}: {cs_line}!r")

        return [
            read_name,
            read_len,
            relative_strand,
            path_start,
            path_end,
            nodes_list,
            orientation_list,
            cs_line,
        ]

    print(f"ERROR: {len(line_elements)} =? {EXPECTED_GAF_TAGS} _ {int(
        line_elements[MAP_Q_ID])} =? {EXPECTED_MAP_Q}",flush=True, file=stderr)
    return None


def parse_cs_tag(cs_string: str):
    """
    This generator iterates over the cs tag string and returns a list of 'steps' that spell the alignment
    between the read and the path.
    It parses the cs tag to represent it as a list of steps. Each step is a tuple representing the
    difference operator (+,-,:,=,*) and the length of the operation in basepairs.
    For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
    I made it as generator at the beginning thinking about walking on the cs on the fly. Probably not useful now.

    Parameters
    ----------
    cs_string : string
        a string spelling the cs tag in the gaf

    Yields
    -------
    The tuple associated to the last cs tag operator.

    """

    # flag characters used to represent the basepair alignment
    # = : identical sequence, spelled [ACGTN]+
    # + : insertion to the reference, spelled [ACGTN]+
    # - : deletion to the reference, spelled [ACGTN]+
    # : : identical sequence, length [0-9]+
    # * : substitution (reference to query) [acgtn][acgtn]

    # example: cs:Z::6724+T:581+A:1027+G:2962-A:278
    # cs:Z::22+T:7+G:21+T:448+T:676-A:666-A:821-A:819-T:455*CA:340+C:192-A:1757+C:947-C:660-C:616+T:315-G:51+T:20*AT:600+T:721+G:82-G:689+T:930-G:88-T:193+GTC:353-G:163-A:297+C:39+C:173-T:368+T:511-T:198-A:615-G:210-C:501-T:648-T:330-C:1211-C:617-GT:9-T:41+T:19-A:244-T:299*GT:5+GCTTT:60+T:172+G:25+A:48+A:22

    flag_chars = ":*+-="
    # 'i' is defined as 5, as the first part of the field is "cs:Z:"
    i = MIN_CS_LEN - 1
    #until 'i' gets to the end of the string
    while i < len(cs_string):
        # if one of the flag that is followed by a sequence string
        if cs_string[i] in "+-=":
            # the flag is the element pointed by 'i' 
            flag_c = cs_string[i]
            count = i
            i += 1
            # scan until you find another tag or the end
            while i < len(cs_string) and cs_string[i] not in flag_chars:
                i += 1
            # report the length of the string (and the flag)
            yield (flag_c, i - count - 1)
            continue
        
        # if it is identical and the next characters define the length of the identity
        elif cs_string[i] == ":":
            number = ""
            i += 1
            # scan until you find another tag or consume the whole string
            while i < len(cs_string) and (cs_string[i] not in flag_chars):
                if i < len(cs_string):
                    number += cs_string[i]
                    i += 1
            # convert to integer the string of the identity length
            yield (":", int(number))
            continue

        # if it is a sustitution, it has length 2
        elif cs_string[i] == "*":
            i = i + 3
            yield ("*", 1)
            continue

        else:
            i = i + 1
