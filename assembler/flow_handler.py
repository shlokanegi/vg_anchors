from assembler.gaf_reader import GafReader
from assembler.anchor_dictionary_builder import SnarlAnchor
from assembler.alignment_processor import AlignAnchor
import assembler.line_parser as lp

class Orchestrator:
    
    def __init__(self, dictionary: dict, graph_path: str, gaf_path: str):
        self.alignment_processor = AlignAnchor(graph_path, dictionary)
        self.gaf_reader = GafReader(gaf_path)

    def process(self):
        print("in process")
        for line in self.gaf_reader.get_lines():
            print("processing line")
            parsed_data = lp.process_line(line)
            print("processed line")
            if parsed_data:
                result = self.alignment_processor.process_alignment(parsed_data)
            # Do something with the result (e.g., print or store)
