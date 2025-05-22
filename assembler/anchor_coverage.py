from collections import defaultdict
from typing import Dict, List, Tuple
import json

class AnchorCoverage:
    def __init__(self):
        self.initial_coverage: Dict[str, int] = defaultdict(int)  # anchor_id -> read count
        self.final_coverage: Dict[str, int] = defaultdict(int)    # anchor_id -> read count
        self.anchor_reads: Dict[str, List[str]] = defaultdict(list)  # anchor_id -> list of read IDs
        
    def record_initial_coverage(self, anchor_id: str, read_id: str):
        """Record initial read coverage for an anchor"""
        self.initial_coverage[anchor_id] += 1
        self.anchor_reads[anchor_id].append(read_id)
        
    def record_final_coverage(self, anchor_id: str, read_id: str):
        """Record final read coverage for an anchor after extension/merging"""
        self.final_coverage[anchor_id] += 1
        
    def get_coverage_stats(self) -> Tuple[Dict[str, int], Dict[str, int]]:
        """Get both initial and final coverage statistics"""
        return self.initial_coverage, self.final_coverage
    
    def save_coverage_stats(self, output_file: str):
        """Save coverage statistics to a JSON file"""
        stats = {
            'initial_coverage': dict(self.initial_coverage),
            'final_coverage': dict(self.final_coverage),
            'anchor_reads': dict(self.anchor_reads)
        }
        with open(output_file, 'w') as f:
            json.dump(stats, f, indent=2)
            
    def load_coverage_stats(self, input_file: str):
        """Load coverage statistics from a JSON file"""
        with open(input_file, 'r') as f:
            stats = json.load(f)
            self.initial_coverage = defaultdict(int, stats['initial_coverage'])
            self.final_coverage = defaultdict(int, stats['final_coverage'])
            self.anchor_reads = defaultdict(list, stats['anchor_reads']) 
