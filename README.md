# VG_ANCHORS

## DOWNLOAD
```
git clone --recursive https://github.com/shlokanegi/vg_anchors.git
cd vg_anchors
```

## Build from Source

### INSTALL
From the `vg_anchors` folder, run the following command to install all dependencies and set up the environment:
```
make init
```

This will:
1. Compile the `sdust` dependency.
2. Compile and install the `bdsg` graph library dependency.
3. Install all required Python packages, including `Click`, `matplotlib`, `seaborn`, `pandas`, `numpy`, `flask`, and `biopython`.
4. Install the `assembler` package itself in an editable mode, suitable for development.

**Prerequisites:**
Before running the installation, please ensure you have the following software installed:
- Python 3.6+
- `pip`
- `make`
- A C/C++ compiler (like `gcc`)
- `cmake`
- `git`

### RUN
Now you can use the tool. 
To build a sentinel to anchor dictionary from the graph use: 
```
vg-anchors build --graph path/to/graph.vg --index path/to/index.dist --output-prefix path/to/output/prefix
```

To get the anchors associated to the alignment to the graph use: 
```
vg-anchors get_anchors --dictionary path/to/dictionary.pkl --graph path/to/graph.vg --alignment path/to/alignment.gaf --fasta path/to/reads.fasta --output path/to/output
```

## Build Executable
To create a standalone executable, run the following commands from the project root:
```bash
# Run from root. It will generate a new config.ini with the latest values from constants.py
python3 generate_config.py

# Build the executable
python3 setup.py build_py
```
This will create an executable file in the `dist/` directory (e.g., `dist/vg-anchors-0.1.0`). You can then run this file directly.


## DEVELOPMENT
For development, the package is installed in editable mode. You can modify the code and the changes will be immediately available without reinstalling.

To run tests:
```bash
python -m unittest discover tests/
```

To start the visualization server (BETA: not in use currently!):
```bash
python server.py
```
Then open http://localhost:8000 in your browser.
