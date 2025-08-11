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
1. Compile the `sdust` dependency (C binary for sequence analysis).
2. Build and install the `libbdsg` dependency (Python bindings for graph operations).
3. Install all required Python packages including:
   - Click (CLI framework)
   - matplotlib, seaborn, pandas, numpy (data visualization and analysis)
   - Flask (web server for visualization interface)
4. Install the assembler package in development mode.

**Prerequisites:**
- Python 3.6+
- pip
- make
- gcc compiler
- cmake (for building libbdsg)
- git

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

**Troubleshooting Installation Conflicts:**
If you encounter errors like "Egg-link does not match installed location" during installation, it usually means there are conflicting installations from different directories. To resolve this:

1. **Check for conflicting installations:**
   ```bash
   pip list | grep -i assembler
   find ~/.local/lib/python*/site-packages -name "*assembler*" -o -name "*vg_assembly*" 2>/dev/null
   ```

2. **Remove conflicting installations:**
   ```bash
   pip uninstall assembler -y
   rm -rf ~/.local/lib/python*/site-packages/assembler*
   ```

3. **Clean any egg-link files pointing to wrong directories:**
   ```bash
   find ~/.local/lib/python*/site-packages -name "*.egg-link" -exec grep -l "vg_assembly" {} \;
   # Remove any egg-link files that point to wrong directories
   ```

4. **Re-run the installation:**
   ```bash
   make init
   ```

**Troubleshooting Module Import Errors:**
If you encounter errors like "dynamic module does not define module export function (PyInit_bdsg)" or similar import errors:

1. **Clean and rebuild the bdsg module:**
   ```bash
   pip uninstall bdsg -y
   cd libbdsg
   rm -rf build/ dist/ *.egg-info/
   pip install -e .
   cd ..
   ```

2. **Or use the clean target:**
   ```bash
   make clean
   make init
   ```

This usually happens when there's a mismatch between the Python version that built the module and the one trying to use it, or when build artifacts are corrupted.

## RUN
Now you can use the tool. 
To build a sentinel to anchor dictionary from the graph use: 
```
vg_anchor build --graph path/to/graph.vg --index path/to/index.dist --output-prefix path/to/output/prefix
```

To get the anchors associated to the alignment to the graph use: 
```
vg_anchor get_anchors --dictionary path/to/dictionary.pkl --graph path/to/graph.vg --alignment path/to/alignment.gaf --fasta path/to/reads.fasta --output path/to/output
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
