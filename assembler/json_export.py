import json

class JsonExporter:
    
    def __init__(self, json_file_path: str):
        self.out_file_path = json_file_path

    def write_lines(self, anchors):
        with open(self.out_file_path, 'w', encoding='utf-8') as f:
            for anchor in anchors:
                print(json.dumps(anchors, f, ensure_ascii=False)+ '\n', file=f)