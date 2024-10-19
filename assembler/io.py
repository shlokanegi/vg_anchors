import sys
import numpy as np

def process_input_gaf():
    """
    Processing the output of a gaf alignment, extracting the path
    with comma-separated node_ids and node_orientations ('+' or '-').

    Returns:
        Tuple:
        - numpy.ndarray: An array containing the node_ids
        - numpy.ndarray: A boolean array, where '+' is True and '-' is false.

    Raises:
        ValueError: If an entry in the data cannot be parsed as an integer or valid sign.
    
    Example:
        Input Line: "1,+,2,-,3,-\n"
        Output:    (array([1, 2, 3]), array([True, False, False])) 
    """
    ids_list = []  # for
    orientations_list = []     #

    for line in sys.stdin:
        # Stipping withespaces/newlines and then splitting the line by '\t'
        fields = line.strip().split('\t')
        
        # if mapQ is not high, do not ingest

        # The sequence name is in the field 0
        sequence_id = fields[0]
        # The path is in the field 5
        path_field = fields[5].split(',')

        for item in path_field:
            item = item.strip()  # Ensure no leading or trailing spaces

            # Instead of storing char for the node orientation, 
            # I store boolean values, should be more memory efficient
            if item == '+':
                ids_list.append(True)  # True for '+'
            elif item == '-':
                orientations_list.append(False)  # False for '-'
            else:
                try:
                    node = int(item)  # Convert to integer
                    ids_list.append(node)
                except ValueError:
                    pass  # Handle any unexpected cases (optional)

    
    # Converting Python lists to NumPy arrays, ight be helpful later in the processing
    int_array = np.array(ids_list, dtype=np.int32)  
    sign_array = np.array(orientations_list, dtype=np.bool_)    

    # Output/return the processed data
    return path