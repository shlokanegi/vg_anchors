class GafReader:
    def __init__(self, text_file_path: str):
        self.file_path = text_file_path

    def get_lines(self):
        with open(self.file_path, 'r') as f:
            for line in f:
                yield line.strip()
