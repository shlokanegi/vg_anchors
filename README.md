#### VG_ANCHORS

##### DOWNLOAD
```
git clone --recursive https://github.com/frankandreace/vg_assembly.git
cd vg_assembly
```

##### INSTALL
From the `vg_assembly` folder, run the following command to install all dependencies and set up the environment:
```
make init
```

This will:
1. Compile the `libbdsg` and `sdust` dependencies.
2. Install the required Python packages.
3. Set up the necessary environment variables.


##### RUN
Now you can use the tool. 
To build a sentinel to anchor dictionary from the graph use: 
```
vg_anchor build --graph path/to/graph.vg --index path/to/index.dist --output-prefix path/to/output/prefix
```

To get the anchors associated to the alignment to the graph use: 
```
vg_anchor get_anchors --dictionary path/to/dictionary.pkl --graph path/to/graph.vg --alignment path/to/alignment.gaf --fasta path/to/reads.fasta --output path/to/output
```
