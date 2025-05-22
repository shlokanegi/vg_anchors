#!/usr/bin/env python3

import unittest
from assembler.aligner import AlignAnchor
from assembler.anchor import Anchor
from assembler.node import Node
from bdsg.bdsg import PackedGraph
import os
import tempfile

class TestAnchorExtension(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        self.aligner = AlignAnchor()
        
        # Create a simple test graph
        self.graph = PackedGraph()
        
        # Create test nodes
        self.node1 = Node(1, 100, True)  # length 100, forward orientation
        self.node2 = Node(2, 150, True)
        self.node3 = Node(3, 200, True)
        self.node4 = Node(4, 120, True)
        
        # Create test anchors
        self.anchor1 = Anchor([self.node1, self.node2])
        self.anchor2 = Anchor([self.node2, self.node3])
        
        # Add test reads
        self.test_reads = [
            ["read1", 0, 0, 100, 50, 20, 30],  # [read_id, strand, start, end, match_limit, cs_left, cs_right]
            ["read2", 0, 0, 100, 50, 20, 30],
            ["read3", 0, 0, 100, 50, 20, 30]
        ]
        
        # Initialize test data
        self.anchor1.bp_matched_reads = self.test_reads
        self.anchor2.bp_matched_reads = self.test_reads
        
        # Set up snarl dictionaries
        self.aligner.snarl_to_anchors_dictionary = {
            "snarl1": [self.anchor1],
            "snarl2": [self.anchor2]
        }
        
        # Create temporary files for graph and dictionary
        self.temp_dir = tempfile.mkdtemp()
        self.graph_path = os.path.join(self.temp_dir, "test_graph.pg")
        self.dict_path = os.path.join(self.temp_dir, "test_dict.pkl")
        
        # Save test graph
        self.graph.serialize(self.graph_path)
        
    def tearDown(self):
        """Clean up test fixtures after each test method."""
        # Remove temporary files
        if os.path.exists(self.graph_path):
            os.remove(self.graph_path)
        if os.path.exists(self.dict_path):
            os.remove(self.dict_path)
        os.rmdir(self.temp_dir)
    
    def test_extending_anchors_by_merging(self):
        """Test merging of anchors from adjacent snarls."""
        # Set up test data
        snarl_ids = ["snarl1", "snarl2"]
        current_snarl_anchors = [self.anchor1]
        anchors_to_discard = []
        
        # Test merging
        new_anchors, updated_idx = self.aligner._extending_anchors_by_merging(
            snarl_ids, 0, "snarl1", "snarl2",
            current_snarl_anchors, True, anchors_to_discard, True
        )
        
        # Verify results
        self.assertIsNotNone(new_anchors)
        self.assertGreater(len(new_anchors), 0)
        self.assertIsInstance(updated_idx, int)
    
    def test_extending_anchors_to_1_degree_node(self):
        """Test extension of anchors into 1-degree nodes."""
        # Set up test data
        current_snarl_anchors = [self.anchor1]
        anchors_to_discard = []
        snarl_ids = ["snarl1"]
        
        # Test extension
        extended_anchors = self.aligner._extending_anchors_to_1_degree_node(
            self.node3, self.node2, "snarl1",
            current_snarl_anchors, True, anchors_to_discard,
            snarl_ids, 0
        )
        
        # Verify results
        self.assertIsNotNone(extended_anchors)
        self.assertGreaterEqual(len(extended_anchors), len(current_snarl_anchors))
    
    def test_try_extension(self):
        """Test the extension process with different iteration parameters."""
        # Set up test data
        current_snarl_anchors = [self.anchor1]
        anchors_to_discard = []
        per_anchor_max_bps = [50]  # 50 bp extension allowed
        
        # Test extension with no drops
        extended_anchors = self.aligner._try_extension(
            current_snarl_anchors, "snarl1", "snarl2",
            anchors_to_discard, per_anchor_max_bps, True, 0
        )
        
        # Verify results
        self.assertIsNotNone(extended_anchors)
        self.assertGreaterEqual(len(extended_anchors), len(current_snarl_anchors))
    
    def test_extending_snarl_boundaries(self):
        """Test the main snarl boundary extension process."""
        # Set up test data
        current_snarl_anchors = [self.anchor1]
        anchors_to_discard = []
        snarl_ids = ["snarl1"]
        
        # Test extension with different iterations
        for iteration in range(3):
            self.aligner._extending_snarl_boundaries(
                current_snarl_anchors, "snarl1",
                snarl_ids, 0, anchors_to_discard, iteration
            )
            
            # Verify results after each iteration
            self.assertIsNotNone(current_snarl_anchors)
            self.assertGreaterEqual(len(current_snarl_anchors), 0)
    
    def test_coverage_tracking(self):
        """Test that coverage is properly tracked during extension."""
        # Record initial coverage
        for read in self.test_reads:
            self.aligner.anchor_coverage.record_initial_coverage(
                f"{self.anchor1!r}", read[0]
            )
        
        # Verify initial coverage
        initial_coverage = self.aligner.anchor_coverage.initial_coverage[f"{self.anchor1!r}"]
        self.assertEqual(initial_coverage, len(self.test_reads))
        
        # Simulate extension and record final coverage
        for read in self.test_reads:
            self.aligner.anchor_coverage.record_final_coverage(
                f"{self.anchor1!r}", read[0]
            )
        
        # Verify final coverage
        final_coverage = self.aligner.anchor_coverage.final_coverage[f"{self.anchor1!r}"]
        self.assertEqual(final_coverage, len(self.test_reads))

if __name__ == '__main__':
    unittest.main()
