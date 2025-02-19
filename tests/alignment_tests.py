import unittest

from assembler.node import Node
from assembler.anchor import Anchor
from assembler.aligner import verify_path_concordance, verify_sequence_agreement


class TestVerifyPathConcordance(unittest.TestCase):

    def setUp(self):
        # defining nodes
        self.node_A = Node(1, 81, True)
        self.node_B = Node(2, 1, True)
        self.node_C = Node(3, 2, True)
        self.node_D = Node(4, 41, True)
        self.node_E = Node(5, 56, True)
        

    def test_concordant_matching_alignment(self):
        # this test verifies if the function returns correct values for a concordant and matching path

        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_D)

        # setting up alignment values

        alignment_node_id_list = [1,2,4]
        alignment_orientation_list = [True, True, True]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

        #EXPECTED CALCULATION:
        # MATCHING: True
        matching = True
        # START_WALK: 200 - 81 + (81 + 1) // 2 = 200 - 81 + 41 = 160
        # END WALK: 200 + 1 + 41 - (41+1) // 2 = 200 + 1 + 41 - 21 = 221
        start = 160
        end = 221
        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)

    def test_out_of_range_concordant(self):
        # this test verifies that the method correctly returns FALSE as the alignment is out_of_range compared to the anchor
        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_D)

        # setting up alignment values

        alignment_node_id_list = [1,2]
        alignment_orientation_list = [True, True]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

        #EXPECTED CALCULATION:
        # MATCHING: FALSE - THE ALIGNMENT STOPS AFTER NODE 2
        matching = False
        # START AND END ARE SET TO 0 IF MATCHING IS FALSE
        start = 0
        end = 0

        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)



    def test_out_of_range_non_concordant(self):
        # this test verifies if the function returns correct values for a concordant and matching path

        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_D)

        # setting up alignment values

        alignment_node_id_list = [4,2]
        alignment_orientation_list = [False, False]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

        #EXPECTED CALCULATION:
        # MATCHING: FALSE - THE ALIGNMENT STOPS AFTER NODE 2
        matching = False
        # START AND END ARE SET TO 0 IF MATCHING IS FALSE
        start = 0
        end = 0

        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)

    def test_concordant_unmatching_alginment(self):
        # this test verifies if the function returns correct values for a concordant and matching path

        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_E)

        # setting up alignment values

        alignment_node_id_list = [1,2,4]
        alignment_orientation_list = [True, True, True]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

        #EXPECTED CALCULATION:
        # MATCHING: FALSE, Node_E is not Node_D
        matching = False
        # START AND END are 0 if FALSE
        start = 0
        end = 0
        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)

    def test_non_concordant_matching_alignment(self):
        # this test verifies if the function returns correct values for a concordant and matching path

        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_E)

        # setting up alignment values

        alignment_node_id_list = [5,2,1]
        alignment_orientation_list = [False, False, False]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

        #EXPECTED CALCULATION:
        # MATCHING: True
        matching = True
        # START_WALK: 200 - 56 + (56 + 1) // 2 = 200 - 56 + 28 = 172
        # END WALK: 200 + 1 + 81 - (81+1) // 2 = 200 + 1 + 81 - 41 = 241
        start = 172
        end = 241
        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)

    def test_non_concordant_not_matching_alignment(self):
        # this test verifies if the function returns correct values for a concordant and matching path

        # adding nodes to anchor
        anchor = Anchor()
        anchor.add(self.node_A)
        anchor.add(self.node_B)
        anchor.add(self.node_E)

        # setting up alignment values

        alignment_node_id_list = [5,2,3]
        alignment_orientation_list = [False, False, False]

        # position in the alignment
        alignment_position = 1
        #id of the sentinel (node B)
        node_id_sentinel = 2
        # length walked in the alignment.
        walked_length = 200

         #EXPECTED CALCULATION:
        # MATCHING: FALSE, Node_C is not Node_E
        matching = False
        # START AND END are 0 if FALSE
        start = 0
        end = 0
        expected_result = (matching, start, end)

        function_result = verify_path_concordance(alignment_position,node_id_sentinel,alignment_node_id_list,alignment_orientation_list,anchor,walked_length)

        self.assertEqual(function_result, expected_result)


class TestVerifySequenceAgreement(unittest.TestCase):

    def setUp(self):
        # defining nodes
        self.node_A = Node(1, 81, True)
        self.node_B = Node(2, 1, True)
        self.node_C = Node(3, 2, True)
        self.node_D = Node(4, 41, True)
        self.node_E = Node(5, 56, True)
    
    def test_correct_alignment(self):

        anchor_bp_start = 0
        anchor_bp_end = 0
        cigar_tags_walk = []
        start_in_path = 0
        end_in_path = 0
        pass

if __name__ == "__main__":
    unittest.main()