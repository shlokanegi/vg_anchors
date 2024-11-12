### TODO LIST - 8 Nov 2024 update

- [ ] Test on 1Mbp in chr20 9 (Range 150k-1150k in CHM13)
- [x] Download Revio HiFi reads
- [x] Avoid reloading the graph when processing gaf. Do not store handles, store only the info (node_id, isReversed, node_length). This could also be more beneficial for initial dictionary constriction (see comparison).
- [x] Serialize dictionary (not uing handles) so that can be loaded instead of recomputed.
- [ ] Understand how to use paths in GBWT instead of paths in PackedGraph 
- [ ] Store reads associated to alingmen in a different dictionrary
- [ ] Write unit tests for cpp implementation
- [ ] Write cpp implementation
- [x] Find a better way than to copy anchors when traversing (to simplify code for backward traversal)
- [ ] Add a dictionary to verify that it is not processing a read twice. (Use only primary alignments - should be guaranteed by the mapq quality though)