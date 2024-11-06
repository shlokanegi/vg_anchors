import numpy as np

from assembler.constants import EXPECTED_GAF_TAGS, EXPECTED_MAP_Q, MAP_Q_ID, READ_NAME_ID, RELATIVE_STRAND_ID, READ_LEN, PATH_START_ID,PATH_END_ID, PATH_ID, CS_TAG_ID


def process_line(gaf_line: str):
    line_elements = gaf_line.split()
    # map_q = int(line_elements[MAP_Q_ID])

    # print(f"Line has {len(line_elements)} tags, with mapq of {map_q}")
    # First verify that the gaf line contains an usable alignment
    if (len(line_elements) == EXPECTED_GAF_TAGS) and int(line_elements[MAP_Q_ID]) == EXPECTED_MAP_Q:

        #extract needed tags
        read_name = line_elements[READ_NAME_ID]
        read_len = int(line_elements[READ_LEN])
        relative_strand = True if line_elements[RELATIVE_STRAND_ID] == "+" else False
        path_start = int(line_elements[PATH_START_ID])
        path_end = int(line_elements[PATH_END_ID])
        # print("extracted tags")
        # decompose the path into nodes and orientations arrays
        nodes_list = []
        orientation_list = []
        # print(line_elements[PATH_ID])

        # have to rewrite this. need to parse every number
        curr_node_string = ""
        for char in line_elements[PATH_ID]:
            if char in"><":
                orientation_list.append(True if char == ">" else False)
                if len(curr_node_string) != 0:
                    nodes_list.append(int(curr_node_string))
                    curr_node_string = ""
            else:
                curr_node_string += char
        if len(curr_node_string) != 0:
            nodes_list.append(int(curr_node_string))
        
        # print("extracted path")
        # decompose the cs tag into alignment steps
        # print("Producing CS line...", end=" ", flush=True)
        if len(line_elements[CS_TAG_ID]) > 6: 
            cs_line = [ i for i in parse_cs_line(line_elements[CS_TAG_ID])]
        else:
            # print(f"cs line not valid, length is {len(line_elements[CS_TAG_ID])}")
            return None
        print(f"{read_name}: {cs_line}!r")
        # print("Done")
        return [read_name, read_len, relative_strand, path_start, path_end, np.array(nodes_list, dtype=np.int32), np.array(orientation_list, dtype=np.bool_), cs_line]
        # return {'read_name': read_name, 'relative_strand' : relative_strand, 'path_start': path_start,\
        #     'path_end' : path_end, 'nodes':np.array(nodes_list, dtype=np.int32), \
        #         'orientations' : np.array(orientation_list, dtype=np.bool_), 'cs' : cs_line}
    return None


def parse_cs_line(cs_string: str):
    """
    This function iterates over the cs tag string and returns a list of 'steps' that spell the alignment 
    between the read and the path.
    It parses the cs tag to represent it as a list of steps. Each step is a tuple representing the 
    difference operator (+,-,:,=,*) and the length of the operation in basepairs.
    For cs tag description see : https://lh3.github.io/minimap2/minimap2.html#10
    """

    flag_chars = ":*+-="
    i = 5
    # print(f"Size: {len(cs_string)}", flush=True)
    # print(cs_string)
    while(i < len(cs_string)):
        # print()
        if cs_string[i] in "+-=":
            flag_c = cs_string[i]
            count = i
            i += 1
            while(i < len(cs_string) and cs_string[i] not in flag_chars):
                i += 1
            # print(f"{flag_c},{i - count - 1}; i:{i}", flush=True)
            yield (flag_c,i - count - 1)
            continue
        
        elif cs_string[i] == ":":
            number = ""
            i += 1
            while(i < len(cs_string) and (cs_string[i] not in flag_chars)):
                if i < len(cs_string):
                    number += cs_string[i]
                    i += 1
            # print(f":,{int(number)}; i:{i}", flush=True)
            yield (":", int(number))
            continue

        elif cs_string[i] == "*":
            i = i + 3
            # print(f"*,1; i:{i}")
            yield ("*",1)
            continue
        
        else: 
            i = i+1