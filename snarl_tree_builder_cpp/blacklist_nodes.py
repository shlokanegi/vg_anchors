# blacklist_nodes.py
# Pipeline to identify and blacklist centromere nodes in the graph
#
# 1. Read the CHM13 censat BED file and identify boundary nodes defining big enough centromere regions
#    STATUS: ✓ COMPLETED
#    - Script: process_censat_regions.py
#    - Input: chm13v2.0.cenSat.v2.0.bed (2799 regions)
#    - Output: censat_merged_regions.bed (324 merged regions)
#    - Gap analysis: censat_gaps_report.txt
#    - Merge threshold: 50 kb (merges adjacent/closely-spaced regions)
#    - Result: Identified 25 large centromere blocks (>1 Mb), ranging from 1-44 Mb
#
# 2. (SKIPPED) Using CHM13 coordinates, find the corresponding nodes in the graph
#    STATUS: SKIPPED (not needed - gbz-query accepts coordinates directly)
#    - Original plan was to use bdsg to convert coordinates to node IDs
#    - Now using gbz-query --interval which accepts coordinate ranges directly
#
# 3. (SKIPPED) Make sure both boundary nodes of a censat window/region are in the same chain.
#    STATUS: SKIPPED (not needed - gbz-query handles this internally)
#
# 4. Extract subgraphs using gbz-query and save all nodes as blacklist.
#    STATUS: ✓ COMPLETED
#    - Script: step4_extract_blacklist_nodes.py
#    - Process:
#      a) Read censat_merged_regions.bed
#      b) For each region, run gbz query with --sample CHM13 --contig <chr> --interval <start>..<end> --snarls <gbz_db>
#      c) Parse GFA output, extract node IDs from 'S' lines
#      d) Combine nodes per chromosome (multiple regions per chromosome possible)
#      e) Save unique, sorted node IDs as binary int64_t files (one per chromosome)
#    - Uses multi-threading for parallel processing
#    - Outputs:
#      * blacklist_nodes/<chr>_blacklist_nodes.bin (per-chromosome binary files)
#      * blacklist_nodes/all_blacklist_nodes.bin (combined file)
#      * blacklist_nodes/blacklist_summary.tsv (summary with counts)
#    - Note: chrY is skipped (excluded from processing)
