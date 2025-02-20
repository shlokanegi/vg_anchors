import unittest

from assembler.node import Node
from assembler.anchor import Anchor
from assembler.aligner import verify_path_concordance, verify_sequence_agreement
from assembler.parser import parse_cs_tag


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
    
    def test_correct_alignment_1(self):

        anchor_bp_start = 0
        anchor_bp_end = 0
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
    
        expected_result = (True,0,0)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)
    
    
    def test_correct_alignment_2(self):

        anchor_bp_start = 0
        anchor_bp_end = -1
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
        
        expected_result = (False,0,0)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_3(self):

        anchor_bp_start = 0
        anchor_bp_end = 1
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
        
        expected_result = (True,0,1)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_4(self):

        anchor_bp_start = 0
        anchor_bp_end = 100
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
        
        expected_result = (True,0,100)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_5(self):

        anchor_bp_start = 99
        anchor_bp_end = 100
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
        
        expected_result = (True,99,100)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_5(self):

        anchor_bp_start = 99
        anchor_bp_end = 101
        cigar_tags_walk = [(":",100)]
        start_in_path = 0
        end_in_path = 100
        
        expected_result = (False,0,0)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_6(self):

        anchor_bp_start = 110
        anchor_bp_end = 120
        cigar_tags_walk = [("+",100),(":",100),("+",10),("-",10),(":",100)]
        start_in_path = 0
        end_in_path = 300
        
        expected_result = (True,210,220)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_7(self):

        anchor_bp_start = 80
        anchor_bp_end = 140
        cigar_tags_walk = [("+",100),(":",100),("-",10),(":",100)]
        start_in_path = 0
        end_in_path = 300
        
        expected_result = (False,0,0)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)

    def test_correct_alignment_8(self):

        anchor_bp_start = 111
        anchor_bp_end = 151
        cigar_tags_walk = [("+",100),(":",100),("-",10),(":",100)]
        start_in_path = 0
        end_in_path = 300
        
        expected_result = (True,201,241)
        function_result = verify_sequence_agreement(anchor_bp_start, anchor_bp_end, cigar_tags_walk, start_in_path, end_in_path)

        self.assertEqual(function_result,expected_result)


class TestParseCsTag(unittest.TestCase):

    def setUp(self):
        return super().setUp()
    
    def test_parsing_example_1(self):
        cigar_string = "cs:Z::6724+T:581+A:1027+G:2962-A:278"

        expected_result = [(":",6724),("+",1),(":",581),("+",1),(":",1027),("+",1),(":",2962),("-",1),(":",278)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_2(self):
        cigar_string = "cs:Z::22+T:7+G:21+T:448+T:676-A:666-A:821-A:819-T:455*CA:340+C:192-A:1757+C:947-C:660-C:616+T:315-G:51+T:20*AT:600+T:721+G:82-G:689+T:930-G:88-T:193+GTC:353-G:163-A:297+C:39+C:173-T:368+T:511-T:198-A:615-G:210-C:501-T:648-T:330-C:1211-C:617-GT:9-T:41+T:19-A:244-T:299*GT:5+GCTTT:60+T:172+G:25+A:48+A:22"

        expected_result = [(":",22),("+",1),(":",7),("+",1),(":",21),("+",1),(":",448),("+",1),(":",676),("-",1),(":",666),("-",1),(":",821),("-",1),(":",819),("-",1),(":",455),("*",1),(":",340),("+",1),(":",192),("-",1),(":",1757),("+",1),(":",947),("-",1),(":",660),("-",1),(":",616),("+",1),(":",315),("-",1),(":",51),("+",1),(":",20),("*",1),(":",600),("+",1),(":",721),("+",1),(":",82),("-",1),(":",689),("+",1),(":",930),("-",1),(":",88),("-",1),(":",193),("+",3),(":",353),("-",1),(":",163),("-",1),(":",297),("+",1),(":",39),("+",1),(":",173),("-",1),(":",368),("+",1),(":",511),("-",1),(":",198),("-",1),(":",615),("-",1),(":",210),("-",1),(":",501),("-",1),(":",648),("-",1),(":",330),("-",1),(":",1211),("-",1),(":",617),("-",2),(":",9),("-",1),(":",41),("+",1),(":",19),("-",1),(":",244),("-",1),(":",299),("*",1),(":",5),("+",5),(":",60),("+",1),(":",172),("+",1),(":",25),("+",1),(":",48),("+",1),(":",22)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_3(self):
        cigar_string = "cs:Z:+AAAA:21:21+AAA-NNNNN:21"

        expected_result = [("+",4),(":",21),(":",21),("+",3),("-",5),(":",21)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_4(self):
        cigar_string = "cs:Z:-AAAA=AA:21+AAA-NNNNN+AA:21-AA"

        expected_result = [("-",4),("=",2),(":",21),("+",3),("-",5),("+",2),(":",21),("-",2)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_5(self):
        cigar_string = "cs:Z:-AAAA"

        expected_result = [("-",4)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)
    
    def test_parsing_example_6(self):
        cigar_string = "cs:Z:+ANNANNNTGATAAA"

        expected_result = [("+",14)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_7(self):
        cigar_string = "cs:Z:=TAGATGATAGATNNCT"

        expected_result = [("=",16)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_8(self):
        cigar_string = "cs:Z:*AT"

        expected_result = [("*",1)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

    def test_parsing_example_9(self):
        cigar_string = "cs:Z::43452"

        expected_result = [(":",43452)]

        function_result = [i for i in parse_cs_tag(cigar_string)]

        self.assertEqual(function_result,expected_result)

if __name__ == "__main__":
    unittest.main()