from assembler.gaf_reader import GafReader
from assembler.anchor_dictionary_builder import SnarlAnchor
from assembler.alignment_processor import AlignAnchor
import assembler.line_parser as lp
import time

class Orchestrator:
    
    def __init__(self, dictionary: dict, graph_path: str, gaf_path: str):
        self.alignment_processor = AlignAnchor(graph_path, dictionary)
        self.gaf_reader = GafReader(gaf_path)

    def process(self):
        # print("in process")
        for line in self.gaf_reader.get_lines():
            # print("processing line")
            t0 = time.time()
            parsed_data = lp.process_line(line)
            t1 = time.time()
            # print("processed line")
            if parsed_data:
                print(f"PROCESSING READ {parsed_data[0]}...", end="", flush=True)
                result = self.alignment_processor.process_alignment(parsed_data)
                print(f"Done in {time.time()-t1}. Parsed in {t1-t0}.")

            # Do something with the result (e.g., print or store)

    def dump_anchors(self, out_file: str):
        self.alignment_processor.dump_valid_anchors(out_file)