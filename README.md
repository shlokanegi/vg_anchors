#### VG_ANCHORS

##### DOWNLOAD
```
git clone --recursive https://github.com/frankandreace/vg_assembly.git
```

##### INSTALL
From the vg_assembly folder do:
```
pip install -e .
```

##### RUN
Now you can use the tool. 
To build a sentinel to anchor dictionary from the graph use: 
```
vg_anchor build --graph path/to/graph.vg --index path/to/index.dist --output-prefix path/to/output/prefix
```

To get the anchors associated to the alignment to the graph use: 
```
vg_anchor get_anchors --dictionary path/to/dictionary.pkl --graph path/to/graph.vg --alignment path/to/alignment.gaf --output path/to/output
```

To verify that the anchors are correct use: 
```
vg_anchor verify-output --anchors path/to/output/anchors.json --fastq reads/used/for/alignment.fastq
```