import unittest
from collections import defaultdict
import sys

from assembler.gaf_reader import GafReader
from assembler.aligner import AlignAnchor
from assembler.builder import AnchorDictionary
import assembler.parser as lp
from assembler.helpers import open_fastq, fastq_lines, fastq_entries, reverse_complement

READ_NAME_POS = 0
ORIENTATION_POS = 1
START_POS = 2
END_POS = 3

class TestVerifyOutputIsCoherent(unittest.TestCase):
    def setUp(self):
        # defining filenames
        self.graph_filename = "data/ont_LC2024_testset.vg"
        self.distance_index_filename = "data/ont_LC2024_testset.dist"
        self.fastq_filename = ["data/ont_LC2024_testset.fastq"]
        self.gaf_filename = "data/ont_LC2024_testset.gaf"

        self.dictionary_builder = AnchorDictionary()
        self.alignment_processor = AlignAnchor()
        self.reads_range_dictionary = defaultdict(list)
        
        #producing dictionary
        self.dictionary_builder.build(self.graph_filename, self.distance_index_filename)
        self.dictionary_builder.fill_anchor_dictionary()
        
        #processing alignment
        self.alignment_processor.ingest(self.dictionary_builder.sentinel_to_anchor, self.graph_filename)
        self.gaf_reader = GafReader(self.gaf_filename)
        
        for _, line in enumerate(self.gaf_reader.get_lines()):
            parsed_data = lp.processGafLine(line)
            self.alignment_processor.processGafLine(parsed_data)
        
        
        
    # def test_anchors_do_not_overlap(self):

    #     for _, anchor_list in self.alignment_processor.anchor_reads_dict.items():
    #         for anchor in anchor_list:
    #             # print(repr(anchor))
    #             for match in anchor:
    #                 sequence_name = match[READ_NAME_POS]
    #                 range_start = int(match[START_POS])
    #                 range_end = int(match[END_POS])
    #                 self.reads_range_dictionary[sequence_name].append(
    #                     [range_start, range_end]
    #                 )
    #     print(f'Reads_range_dict contains {len(self.reads_range_dictionary)} elements')

    #     # This test verifies anchors do not overlap inside the reads.
    #     assertion_counter: int = 0
    #     for _, values in self.reads_range_dictionary.items():
    #         anchors_ranges = []
    #         for start, end in values:
    #             anchors_ranges.append((start,end))
    #         sorted_ranges = sorted(anchors_ranges, key=lambda y: y[0])
    #         curr_end = -1

    #         for start, end in sorted_ranges:
    #             self.assertGreaterEqual(start,curr_end)
    #             assertion_counter += 1
    #             curr_end = end
    #         print(f"Tot assertion computed: {assertion_counter}",flush=True)


    def test_anchors_sequences_match(self):
        # This test verifies anchors do not overlap inside the reads.

        reads_dict = dict()
        for seq_dict in fastq_entries(fastq_lines(self.fastq_filename)):
            print(seq_dict['header'])
            reads_dict[seq_dict['header']] = seq_dict['sequence']
        # for el in reads_dict:
        #     print(el)

        for _, anchor_list in self.alignment_processor.anchor_reads_dict.items():
            to_check_list = []
            for anchor in anchor_list:
                if len(anchor) > 1:
                    for match in anchor:
                        print(repr(match))
                        print(repr(match[READ_NAME_POS]))
                        print(repr(reads_dict[match[READ_NAME_POS]]))
                        sys.exit(1)
                        sequence = reads_dict[match[READ_NAME_POS]]['sequence']
                        reversed_sequence = reverse_complement(sequence)
                        if match[ORIENTATION_POS] == 1:
                            to_check_list.append(reversed_sequence[match[START_POS]:match[END_POS]])
                        else:
                            to_check_list.append(sequence[match[START_POS]:match[END_POS]])
                    anchor_sequence = to_check_list[0]
                    for sequence in to_check_list[1:]:
                        self.assertEqual(len(anchor_sequence),len(sequence))


if __name__ == "__main__":
    unittest.main()