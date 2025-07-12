document.addEventListener('DOMContentLoaded', function() {

    // const infoContent = document.getElementById('info-content');
    const readIdInput = document.getElementById('read-id-input');
    const visualizeReadBtn = document.getElementById('visualize-read-btn');
    const resetViewBtn = document.getElementById('reset-view-btn');
    const increaseNodeSizeBtn = document.getElementById('increase-node-size');
    const decreaseNodeSizeBtn = document.getElementById('decrease-node-size');
    const increaseEdgeWidthBtn = document.getElementById('increase-edge-width');
    const decreaseEdgeWidthBtn = document.getElementById('decrease-edge-width');
    const anchorIdInput = document.getElementById('anchor-id-input');
    const viewAnchorBtn = document.getElementById('view-anchor-btn');

    let nodeSize = 8;
    let edgeWidthMultiplier = 1.0;

    // Load all data files
    Promise.all([
        fetch('snarl_to_anchors_dictionary.json').then(res => res.json()),
        fetch('anchor_to_anchor_connectivity_info.json').then(res => res.json()),
        fetch('anchors_read_info.json').then(res => res.json())
    ]).then(([snarls, connectivity, anchorReads]) => {
        
        console.log("Data loaded successfully");
        console.log("Snarls:", snarls);
        console.log("Connectivity:", connectivity);
        console.log("Anchor Reads:", anchorReads);

        // --- 1. Process Data ---

        const readToAnchors = {};
        for (const anchorId in anchorReads) {
            anchorReads[anchorId].forEach(read => {
                if (!readToAnchors[read.read_id]) {
                    readToAnchors[read.read_id] = [];
                }
                readToAnchors[read.read_id].push(anchorId);
            });
        }

        const nodes = [];
        const snarlIds = Object.keys(snarls).sort((a, b) => a - b);
        const x_spacing = 200;
        const y_spacing = 100;

        snarlIds.forEach((snarlId, snarlIndex) => {
            const anchorsInSnarl = snarls[snarlId];
            anchorsInSnarl.forEach((anchorId, anchorIndex) => {
                nodes.push({
                    group: 'nodes',
                    data: { id: anchorId, snarl: snarlId },
                    position: {
                        x: snarlIndex * x_spacing,
                        y: anchorIndex * y_spacing
                    }
                });
            });
        });

        console.log("Processed Nodes:", nodes);

        const edges = [];
        for (const anchor1 in connectivity) {
            for (const anchor2 in connectivity[anchor1]) {
                edges.push({
                    group: 'edges',
                    data: {
                        id: `${anchor1}-${anchor2}`,
                        source: anchor1,
                        target: anchor2,
                        weight: connectivity[anchor1][anchor2].common_read_count
                    }
                });
            }
        }
        
        console.log("Processed Edges:", edges);

        // --- 2. Initialize Cytoscape ---

        const cy = cytoscape({
            container: document.getElementById('cy'),
            elements: [...nodes, ...edges],
            layout: {
                name: 'preset'
            },
            style: [
                {
                    selector: 'node',
                    style: {
                        'background-color': '#666',
                        'width': nodeSize,
                        'height': nodeSize,
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width': `mapData(weight, 1, 20, ${1 * edgeWidthMultiplier}, ${10 * edgeWidthMultiplier})`, // Map edge weight to thickness
                        'line-color': '#ccc',
                        'target-arrow-color': '#ccc',
                        'target-arrow-shape': 'triangle',
                        'curve-style': 'bezier'
                    }
                },
                {
                    selector: '.highlighted',
                    style: {
                        'background-color': '#337ab7',
                        'line-color': '#337ab7',
                        'target-arrow-color': '#337ab7',
                        'transition-property': 'background-color, line-color, target-arrow-color',
                        'transition-duration': '0.5s'
                    }
                },
                {
                    selector: '.highlighted-red',
                    style: {
                        'background-color': 'red',
                        'border-color': 'red',
                        'transition-property': 'background-color, border-color',
                        'transition-duration': '0.5s'
                    }
                },
                 {
                    selector: '.faded',
                    style: { 'opacity': 0.25 }
                }
            ],
            zoom: 1,
            pan: { x: 0, y: 0 },
            minZoom: 0.1,
            maxZoom: 5
        });

        // --- 3. Event Listeners ---

        // Tooltip logic
        let tippies = {};

        function makeTippy(node) {
            let ref = node.popperRef();
            let dummy = document.createElement('div');
            
            let tip = tippy(dummy, {
                getReferenceClientRect: ref.getBoundingClientRect,
                trigger: 'manual',
                content: () => {
                    let content = document.createElement('div');
                    content.innerHTML = node.id();
                    return content;
                }
            });

            return tip;
        }

        function setTippy(node) {
            if (tippies[node.id()]) {
                tippies[node.id()].destroy();
            }

            tippies[node.id()] = makeTippy(node);
        }

        cy.nodes().forEach(setTippy);

        // Edge Tooltip logic
        let edgeTippies = {};

        function makeEdgeTippy(edge) {
            let ref = edge.popperRef();
            let dummy = document.createElement('div');
            
            let tip = tippy(dummy, {
                getReferenceClientRect: ref.getBoundingClientRect,
                trigger: 'manual',
                content: () => {
                    let content = document.createElement('div');
                    content.innerHTML = `Common Reads: ${edge.data('weight')}`;
                    return content;
                }
            });

            return tip;
        }

        function setEdgeTippy(edge) {
            if (edgeTippies[edge.id()]) {
                edgeTippies[edge.id()].destroy();
            }
            edgeTippies[edge.id()] = makeEdgeTippy(edge);
        }

        cy.edges().forEach(setEdgeTippy);


        cy.on('mouseover', 'node', (e) => {
            const node = e.target;
            if (tippies[node.id()]) {
                 tippies[node.id()].show();
            }
        });

        cy.on('mouseout', 'node', (e) => {
            const node = e.target;
            if (tippies[node.id()]) {
                 tippies[node.id()].hide();
            }
        });

        cy.on('mouseover', 'edge', (e) => {
            const edge = e.target;
            if (edgeTippies[edge.id()]) {
                 edgeTippies[edge.id()].show();
            }
        });

        cy.on('mouseout', 'edge', (e) => {
            const edge = e.target;
            if (edgeTippies[edge.id()]) {
                 edgeTippies[edge.id()].hide();
            }
        });

        cy.on('tap', 'node', function(evt){
            const node = evt.target;
            const anchorId = node.id();
            const reads = anchorReads[anchorId];
            showNodeInfo(anchorId, reads);
        });

        cy.on('tap', 'edge', function(evt){
            const edge = evt.target;
            const sourceId = edge.source().id();
            const targetId = edge.target().id();
            const sourceReads = new Set(anchorReads[sourceId].map(r => r.read_id));
            const targetReads = anchorReads[targetId];
            const commonReads = targetReads.filter(r => sourceReads.has(r.read_id));
            showEdgeInfo(sourceId, targetId, commonReads);
        });

        visualizeReadBtn.addEventListener('click', () => {
            const readId = readIdInput.value.trim();
            if (readId && readToAnchors[readId]) {
                highlightReadJourney(readId);
            } else {
                alert('Read ID not found!');
            }
        });

        readIdInput.addEventListener('keyup', function(event) {
            if (event.key === 'Enter') {
                visualizeReadBtn.click();
            }
        });

        resetViewBtn.addEventListener('click', () => {
            resetHighlight();
            cy.fit();
            // infoContent.innerHTML = '<p>Click on a node or edge to see details.</p>';
        });
        
        increaseNodeSizeBtn.addEventListener('click', () => {
            nodeSize += 2;
            updateStyles();
        });

        decreaseNodeSizeBtn.addEventListener('click', () => {
            nodeSize = Math.max(2, nodeSize - 2);
            updateStyles();
        });

        increaseEdgeWidthBtn.addEventListener('click', () => {
            edgeWidthMultiplier += 0.2;
            updateStyles();
        });

        decreaseEdgeWidthBtn.addEventListener('click', () => {
            edgeWidthMultiplier = Math.max(0.2, edgeWidthMultiplier - 0.2);
            updateStyles();
        });

        viewAnchorBtn.addEventListener('click', () => {
            const anchorId = anchorIdInput.value.trim();
            if (!anchorId) return;

            const node = cy.getElementById(anchorId);
            if (node.length > 0) {
                cy.animate({
                    fit: {
                        eles: node,
                        padding: 200
                    },
                    duration: 500
                });
                node.addClass('highlighted-red');
                setTimeout(() => {
                    node.removeClass('highlighted-red');
                }, 1500);
            } else {
                alert('Anchor ID not found!');
            }
        });

        anchorIdInput.addEventListener('keyup', function(event) {
            if (event.key === 'Enter') {
                viewAnchorBtn.click();
            }
        });


        // --- 4. Helper Functions ---

        function updateStyles() {
            cy.style()
                .selector('node')
                    .style('width', nodeSize)
                    .style('height', nodeSize)
                .selector('edge')
                    .style('width', `mapData(weight, 1, 20, ${1 * edgeWidthMultiplier}, ${10 * edgeWidthMultiplier})`)
                .update();
        }

        function showNodeInfo(anchorId, reads) {
            let html = `<html><head><title>Anchor: ${anchorId}</title><style>body{font-family:sans-serif} table{width:100%; border-collapse:collapse} th,td{border:1px solid #ddd; padding:8px}</style></head><body>`;
            html += `<h3>Anchor: ${anchorId}</h3><h4>Reads:</h4>`;
            if (reads && reads.length > 0) {
                html += '<table><tr><th>Read ID</th><th>Strand</th><th>Start</th><th>End</th></tr>';
                reads.forEach(r => {
                    html += `<tr><td>${r.read_id}</td><td>${r.read_strand}</td><td>${r.read_start}</td><td>${r.read_end}</td></tr>`;
                });
                html += '</table>';
            } else {
                html += '<p>No read information available.</p>';
            }
            html += '</body></html>';
            const newTab = window.open();
            newTab.document.write(html);
            newTab.document.close();
        }

        function showEdgeInfo(sourceId, targetId, commonReads) {
            let html = `<html><head><title>Edge: ${sourceId} - ${targetId}</title><style>body{font-family:sans-serif} table{width:100%; border-collapse:collapse} th,td{border:1px solid #ddd; padding:8px}</style></head><body>`;
            html += `<h3>Edge: ${sourceId} - ${targetId}</h3><h4>Common Reads:</h4>`;
            if (commonReads && commonReads.length > 0) {
                html += '<table><tr><th>Read ID</th><th>Strand</th><th>Start</th><th>End</th></tr>';
                commonReads.forEach(r => {
                    html += `<tr><td>${r.read_id}</td><td>${r.read_strand}</td><td>${r.read_start}</td><td>${r.read_end}</td></tr>`;
                });
                html += '</table>';
            } else {
                html += '<p>No common reads.</p>';
            }
            html += '</body></html>';
            const newTab = window.open();
            newTab.document.write(html);
            newTab.document.close();
        }
        
        function resetHighlight() {
            cy.elements().removeClass('highlighted').removeClass('faded');
        }

        function highlightReadJourney(readId) {
            resetHighlight();
            
            const anchorsForRead = readToAnchors[readId];
            if (!anchorsForRead) return;

            cy.elements().addClass('faded');

            const nodesToHighlight = cy.nodes().filter(node => anchorsForRead.includes(node.id()));
            nodesToHighlight.removeClass('faded').addClass('highlighted');

            // Highlight edges between the highlighted nodes
            for (let i = 0; i < anchorsForRead.length - 1; i++) {
                for (let j = i + 1; j < anchorsForRead.length; j++) {
                    const source = anchorsForRead[i];
                    const target = anchorsForRead[j];
                    // Check both directions
                    const edge1 = cy.edges(`[source = "${source}"][target = "${target}"]`);
                    const edge2 = cy.edges(`[source = "${target}"][target = "${source}"]`);
                    if(edge1) edge1.removeClass('faded').addClass('highlighted');
                    if(edge2) edge2.removeClass('faded').addClass('highlighted');
                }
            }
        }
        
        cy.fit(); // Fit graph to view on load

    }).catch(error => {
        console.error('Error loading data:', error);
        // infoContent.innerHTML = '<p style="color: red;">Error loading visualization data. Please check the console.</p>';
    });
});
