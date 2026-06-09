# VG_ANCHORS

## DOWNLOAD
```
git clone --recursive https://github.com/shlokanegi/vg_anchors.git
cd vg_anchors
```

## INSTALL
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

**Alternative installation:**
If you prefer to install dependencies manually:
```bash
# Install Python dependencies
pip install -r requirements.txt

# Build and install libbdsg
cd libbdsg
pip install -e .
cd ..

# Install the assembler package
pip install -e .
```

## RUN

### Setup config
All constants required for running vg-anchors should be specified in a config file. The default: `config.ini`.

### vg-anchors build
Step 1: Candidate anchor generation from paths in leaf snarls
* Create leaf snarl dictionaries (initial boundaries of anchors)
* Using paths in the graph, generate initial candidate anchors
```
vg-anchors --config /path/to/config.ini build --graph /path/to/graph.pg.vg --index path/to/index.dist --output-prefix path/to/output/prefix
```

### vg-anchors get-anchors
Step 2: Primary anchor building using aligned reads

Step 3: Reliable snarl finding with phasing consistency to neighbouring snarls 
* Identify “linked snarls” based on shared reads
* Snarl compatibility check based on read partitions of shared reads

Step 4: Anchor extension and merging using aligned reads
* Base-level extension within snarl boundary and into next 1-degree node(s)
* Snarl merging
* Independent anchor extension
```
vg_anchor --config /path/to/config.ini get_anchors --dictionary path/to/dictionary.pkl --threads 8 --graph path/to/graph.pg.vg --alignment path/to/alignment.gaf --fasta path/to/reads.fasta --output path/to/output_json
```

## DEVELOPMENT
For development, the package is installed in editable mode. You can modify the code and the changes will be immediately available without reinstalling.

To run tests:
```bash
python -m unittest discover tests/
```

To start the visualization server:
```bash
python server.py
```
Then open http://localhost:8000 in your browser.

## Build a standalone executable
To produce a single-file executable with all dependencies (including `sdust`, `bdsg`, and shared libraries) bundled:

```bash
bash scripts/build_executable.sh
```

The executable will be at `dist/vg-anchors-0.1.0`. You can run it directly:

```bash
./dist/vg-anchors-0.1.0 --help
```

Notes:
- Build runs in an isolated virtual environment under `.venv` to ensure reproducibility.
- `sdust` is compiled from `third_party/sdust` into `bin/sdust`.
- `libbdsg` is built with CMake and installed into the venv; its `.so` files are bundled under `lib/` in the executable.
- The PyInstaller spec `vg-anchors-0.1.0.spec` explicitly collects `bdsg` and sets a runtime hook (`pyi_rth_vg_anchors_libpath.py`) to add the bundled `lib/` to `LD_LIBRARY_PATH` at runtime.


## DEBUG with VS-Code on UCSC Phoenix cluster
On the local terminal, run ssh phoenix-23-debug (currently, only the phoenix-23 node has been set in `~/.ssh/config` for automatic tunneling and port forwarding. However, you can set another compute node in the same way in the `~/.ssh/config` file.)

TODO: Add more notes


## PROFILE

### Line Profiling (Time)
The codebase is already configured with `@profile` decorators on key functions. The `line_profiler` package is included in `setup.py`, so it should be installed with the package. If not, install it manually:
```bash
pip install line_profiler
```

**Usage:**
To profile a vg-anchors command, use `kernprof` (provided by `line_profiler`) with the profiling wrapper script:
```bash
# Profile the 'build' command
kernprof -l scripts/profile_vg_anchors.py --config /path/to/config.ini build \
    --graph /path/to/graph.pg.vg \
    --index /path/to/index.dist \
    --output-prefix /path/to/output/prefix

# Profile the 'get-anchors' command
kernprof -l scripts/profile_vg_anchors.py --config /path/to/config.ini get-anchors \
    --dictionary /path/to/dictionary.pkl \
    --graph /path/to/graph.pg.vg \
    --alignment /path/to/alignment.gaf \
    --fasta /path/to/reads.fasta \
    --output /path/to/output_json \
    --threads 8

# Profile the 'chunker-build' command
kernprof -l scripts/profile_vg_anchors.py --config /path/to/config.ini chunker-build \
    --graph /path/to/graph.pg.vg \
    --index /path/to/index.dist \
    --output-prefix /path/to/output/prefix
```

After running the profiled command, a `.lprof` file will be generated in the current directory (typically `profile_vg_anchors.py.lprof` or similar). To view the results:

```bash
python -m line_profiler profile_vg_anchors.py.lprof > time_profile_results.txt
```

### Function Profiling (cProfile)

Function-level profiling with `cProfile` will be documented here in a future update. This will provide function-level timing statistics.
