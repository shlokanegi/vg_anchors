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
Now you can use it with 
```
vg_anchor build --graph path/to/graph.vg --index path/to/index.dist --output-prefix path/to/output/prefix
vg_anchor get_anchors --dictionary path/to/dictionary.pkl --graph path/to/graph.vg --alignment path/to/alignment.gaf --output path/to/output 
```