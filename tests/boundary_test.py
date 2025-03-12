import unittest
from sys import stderr
from assembler.builder import AnchorDictionary
from assembler.constants import REVERSE_DICTIONARY, FORWARD_DICTIONARY

FILES_PATH = "data/boundary_tests/"



class TestVerifyBoundaryExtension(unittest.TestCase):

    def verify_dictionary_concordance(self, expected, computed):
        self.assertEqual(len(expected),len(computed))
        for key, expected_end in expected.items():
            computed_end = computed[key][0]
            self.assertEqual(expected_end, computed_end)

    def setUp(self):
        self.dictionary_builder = AnchorDictionary()
        

    def test_extension_left(self):
        test_filename = "test_go_left"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder = AnchorDictionary()
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward = {1:5, 5:8}
        expected_boundary_backward = {5:1, 8:5}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )

    def test_extension_right(self):
        test_filename = "test_go_right"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward  = {1:4, 4:8, 8:11}
        expected_boundary_backward = {4:1, 8:4, 11:8}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )


    def test_stop_on_degree_greater_1(self):
        test_filename = "test_stops_on_deg_gt_1"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward  = {1:4, 4:7, 7:10}
        expected_boundary_backward = {4:1, 7:4, 10:7}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )


    def test_goes_correct_left(self):
        test_filename = "test_goes_correct_left"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward  = {1:4, 4:8, 8:12}
        expected_boundary_backward = {4:1, 8:4, 12:8}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )

    def test_goes_correct_right(self):
        test_filename = "test_goes_correct_right"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward  = {1:5, 5:9, 9:12}
        expected_boundary_backward = {12:9, 9:5, 5:1}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )

    def test_stops_correct(self):
        test_filename = "test_stops_correct"
        graph_filename =  FILES_PATH + test_filename + ".vg"
        distance_index_filename = FILES_PATH + test_filename + ".dist"
        print(f"\n----\nTESTING {test_filename}\n-----",flush=True)

        # load files
        self.dictionary_builder.build(graph_filename, distance_index_filename)
        
        # process snarls
        self.dictionary_builder.process_snarls()

        # get boundaries
        extend = True
        self.dictionary_builder.fill_anchor_dictionary(extend)
        self.dictionary_builder.print_anchor_boundaries_stderr()

        # verify correctness
        expected_boundary_forward  = {2:6}
        expected_boundary_backward = {6:2}

        self.verify_dictionary_concordance(expected_boundary_forward, self.dictionary_builder.snarl_boundaries[FORWARD_DICTIONARY] )

        self.verify_dictionary_concordance(expected_boundary_backward, self.dictionary_builder.snarl_boundaries[REVERSE_DICTIONARY] )



if __name__ == "__main__":
    unittest.main()