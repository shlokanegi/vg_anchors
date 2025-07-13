document.addEventListener('DOMContentLoaded', function() {

    // --- CONFIGURATION ---
    const useBundledEdges = true; // Set to false to revert to the old style (one thick edge per pair)
    
    // DOM Elements
    const readIdInput = document.getElementById('read-id-input');
    const visualizeReadBtn = document.getElementById('visualize-read-btn');
    const resetViewBtn = document.getElementById('reset-view-btn');
    const increaseNodeSizeBtn = document.getElementById('increase-node-size');
    const decreaseNodeSizeBtn = document.getElementById('decrease-node-size');
    const increaseEdgeWidthBtn = document.getElementById('increase-edge-width');
    const decreaseEdgeWidthBtn = document.getElementById('decrease-edge-width');
    const anchorIdInput = document.getElementById('anchor-id-input');
    const viewAnchorBtn = document.getElementById('view-anchor-btn');
    const anchorDropInput = document.getElementById('anchor-drop-input');
    const dropAnchorBtn = document.getElementById('drop-anchor-btn');
    const readDropInput = document.getElementById('read-drop-input');
    const dropReadBtn = document.getElementById('drop-read-btn');
    const droppedAnchorsListDiv = document.getElementById('dropped-anchors-list');
    const droppedReadsListDiv = document.getElementById('dropped-reads-list');
    const undoBtn = document.getElementById('undo-btn');
    const saveBtn = document.getElementById('save-btn');
    
    // Cytoscape and style variables
    let cy;
    let nodeSize = 8;
    let edgeWidthMultiplier = 1.0;

    // Data state
    let originalSnarls, originalAnchorReads, originalReadInfo;
    let allReads = new Set();
    let allAnchors = new Set();
    let droppedAnchors = new Set();
    let droppedReads = new Set();
    let history = [];


    // Load all data files
    Promise.all([
        fetch('snarl_to_anchors_dictionary.json').then(res => res.json()),
        fetch('anchors_read_info.json').then(res => res.json()),
        fetch('read_info.json').then(res => res.json())
    ]).then(([snarls, anchorReads, readInfo]) => {
        
        // Store original data
        originalSnarls = snarls;
        originalAnchorReads = anchorReads;
        originalReadInfo = readInfo;

        // Collect all possible IDs for validation
        allAnchors = new Set(Object.keys(originalAnchorReads));
        allReads = new Set(Object.keys(originalReadInfo));

        updateDroppedListsUI();
        updateGraph();

    }).catch(error => {
        console.error('Error loading data:', error);
    });

    function updateDroppedListsUI() {
        droppedAnchorsListDiv.innerHTML = '';
        [...droppedAnchors].sort().forEach(id => {
            const item = document.createElement('div');
            item.className = 'restore-item';
            item.innerHTML = `<span>${id}</span><button data-type="anchor" data-id="${id}">Restore</button>`;
            droppedAnchorsListDiv.appendChild(item);
        });

        droppedReadsListDiv.innerHTML = '';
        [...droppedReads].sort().forEach(id => {
            const item = document.createElement('div');
            item.className = 'restore-item';
            item.innerHTML = `<span>${id}</span><button data-type="read" data-id="${id}">Restore</button>`;
            droppedReadsListDiv.appendChild(item);
        });
    }

    function updateGraph() {
        // --- 1. Process Data based on drops ---
        const nodes = [];
        const snarlIds = Object.keys(originalSnarls).sort((a, b) => {
            const idA = parseInt(a.split('-')[0].replace('snarl_id', ''));
            const idB = parseInt(b.split('-')[0].replace('snarl_id', ''));
            return idA - idB;
        });

        snarlIds.forEach((snarlId, snarlIndex) => {
            const anchorsInSnarl = originalSnarls[snarlId];
            anchorsInSnarl.forEach((anchorId, anchorIndex) => {
                // TODO: Where are we returning back to?
                if (droppedAnchors.has(anchorId)) return;
                nodes.push({
                    group: 'nodes',
                    data: { id: anchorId, snarl: snarlId },
                    // TODO: Why are we multiplying by 200 and 100? Where are these x and y coordinates being used?
                    position: {
                        x: snarlIndex * 200,
                        y: anchorIndex * 100
                    }
                });
            });
        });

        const edges = [];
        if (useBundledEdges) {
            // --- New "Bundled" Edge Logic ---
            for (const readId in originalReadInfo) {
                if (droppedReads.has(readId)) continue;
                
                const journey = originalReadInfo[readId].journey;
                for (let i = 0; i < journey.length - 1; i++) {
                    const source = journey[i];
                    const target = journey[i+1];

                    if (droppedAnchors.has(source) || droppedAnchors.has(target)) continue;

                    edges.push({
                        group: 'edges',
                        data: {
                            id: `edge-${readId}-${source}-${target}`, // Unique ID per read connection
                            source: source,
                            target: target,
                            read_id: readId // Store the read ID for click events
                        }
                    });
                }
            }
        } else {
            // --- Old "Weighted" Edge Logic ---
            const edgeWeights = {}; 
            for (const readId in originalReadInfo) {
                if (droppedReads.has(readId)) continue;
                const journey = originalReadInfo[readId].journey;
                for (let i = 0; i < journey.length - 1; i++) {
                    const source = journey[i];
                    const target = journey[i+1];
                    if (droppedAnchors.has(source) || droppedAnchors.has(target)) continue;
                    const key = [source, target].sort().join('-');
                    edgeWeights[key] = (edgeWeights[key] || 0) + 1;
                }
            }
            for (const key in edgeWeights) {
                const [source, target] = key.split('-');
                edges.push({
                    group: 'edges',
                    data: {
                        id: `${source}-${target}`,
                        source: source,
                        target: target,
                        weight: edgeWeights[key]
                    }
                });
            }
        }
        
        // Filter out isolated nodes
        const nodesWithEdges = new Set();
        edges.forEach(edge => {
            nodesWithEdges.add(edge.data.source);
            nodesWithEdges.add(edge.data.target);
        });
        const filteredNodes = nodes.filter(node => nodesWithEdges.has(node.data.id));

        // --- 2. Initialize or Update Cytoscape ---
        if (cy) {
            cy.destroy();
        }

        cy = cytoscape({
            container: document.getElementById('cy'),
            elements: [...filteredNodes, ...edges],
             style: [
                {
                    selector: 'node',
                    style: {
                        'background-color': '#000000',
                        'width': nodeSize,
                        'height': nodeSize,
                    }
                },
                useBundledEdges ? { // Style for "Bundled" mode
                    selector: 'edge',
                    style: {
                        'width': 1.5,
                        'line-color': '#cccccc',
                        'curve-style': 'unbundled-bezier'
                    }
                } : { // Style for "Weighted" mode
                    selector: 'edge',
                    style: {
                        'width': `mapData(weight, 1, 20, ${1 * edgeWidthMultiplier}, ${10 * edgeWidthMultiplier})`,
                        'line-color': '#cccccc',
                        'curve-style': 'bezier',
                        'control-point-distance': '40'
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
            layout: { name: 'preset' },
            zoom: 1,
            pan: { x: 0, y: 0 },
            minZoom: 0.1,
            maxZoom: 5
        });

        // --- 3. Re-add Event Listeners ---
        setupEventListeners();
        cy.fit();
    }


    function setupEventListeners() {
        // Tooltip logic
        let tippies = {};
        function makeTippy(node) {
            let ref = node.popperRef();
            let dummy = document.createElement('div');
            let tip = tippy(dummy, { getReferenceClientRect: ref.getBoundingClientRect, trigger: 'manual', content: () => { let div = document.createElement('div'); div.innerHTML = node.id(); return div; } });
            return tip;
        }
        function setTippy(node) { if (tippies[node.id()]) { tippies[node.id()].destroy(); } tippies[node.id()] = makeTippy(node); }
        cy.nodes().forEach(setTippy);
        cy.on('mouseover', 'node', (e) => tippies[e.target.id()]?.show());
        cy.on('mouseout', 'node', (e) => tippies[e.target.id()]?.hide());

        // Edge Tooltip logic
        let edgeTippies = {};
        function makeEdgeTippy(edge) {
            let ref = edge.popperRef();
            let dummy = document.createElement('div');
            let tip = tippy(dummy, { getReferenceClientRect: ref.getBoundingClientRect, trigger: 'manual', content: () => { let div = document.createElement('div'); div.innerHTML = `Common Reads: ${edge.data('weight')}`; return div; } });
            return tip;
        }
        function setEdgeTippy(edge) { if (edgeTippies[edge.id()]) { edgeTippies[edge.id()].destroy(); } edgeTippies[edge.id()] = makeEdgeTippy(edge); }
        cy.edges().forEach(setEdgeTippy);
        cy.on('mouseover', 'edge', (e) => edgeTippies[e.target.id()]?.show());
        cy.on('mouseout', 'edge', (e) => edgeTippies[e.target.id()]?.hide());


        cy.on('tap', 'node', function(evt){
            const node = evt.target;
            const anchorId = node.id();
            const reads = originalAnchorReads[anchorId].filter(r => !droppedReads.has(r.read_id));
            showNodeInfo(anchorId, reads);
        });

        cy.on('tap', 'edge', function(evt){
            const edge = evt.target;
            const sourceId = edge.source().id();
            const targetId = edge.target().id();

            if (useBundledEdges) {
                const readId = edge.data('read_id');
                const readDetails = originalAnchorReads[sourceId]?.find(r => r.read_id === readId) || originalAnchorReads[targetId]?.find(r => r.read_id === readId);
                showEdgeInfo(sourceId, targetId, readDetails ? [readDetails] : []);
            } else {
                // Re-calculate common reads on-the-fly for the info panel
                const commonReads = [];
                for (const readId in originalReadInfo) {
                    if (droppedReads.has(readId)) continue;
                    const journey = originalReadInfo[readId].journey;
                    for (let i = 0; i < journey.length - 1; i++) {
                         // Check for both directions
                        if ((journey[i] === sourceId && journey[i+1] === targetId) || (journey[i] === targetId && journey[i+1] === sourceId)) {
                            const readInfoForEdge = originalAnchorReads[journey[i]]?.find(r => r.read_id === readId) || originalAnchorReads[journey[i+1]]?.find(r => r.read_id === readId);
                            if (readInfoForEdge) {
                                commonReads.push(readInfoForEdge);
                            }
                        }
                    }
                }
                showEdgeInfo(sourceId, targetId, commonReads);
            }
        });
    }

    // --- Global Control Listeners ---
    visualizeReadBtn.addEventListener('click', () => {
        const readId = readIdInput.value.trim();
        // This logic needs to be adapted for the new data flow, for now it will highlight on the current view
        if (readId && allReads.has(readId) && !droppedReads.has(readId)) {
            highlightReadJourney(readId);
        } else {
            alert('Read ID not found or has been dropped!');
        }
    });

    readIdInput.addEventListener('keyup', (e) => { if (e.key === 'Enter') visualizeReadBtn.click(); });
    resetViewBtn.addEventListener('click', () => { resetHighlight(); cy.fit(); });
    increaseNodeSizeBtn.addEventListener('click', () => { nodeSize += 2; updateStyles(); });
    decreaseNodeSizeBtn.addEventListener('click', () => { nodeSize = Math.max(2, nodeSize - 2); updateStyles(); });
    increaseEdgeWidthBtn.addEventListener('click', () => { edgeWidthMultiplier += 0.2; updateStyles(); });
    decreaseEdgeWidthBtn.addEventListener('click', () => { edgeWidthMultiplier = Math.max(0.2, edgeWidthMultiplier - 0.2); updateStyles(); });
    
    viewAnchorBtn.addEventListener('click', () => {
        const anchorId = anchorIdInput.value.trim();
        if (!anchorId) return;
        const node = cy.getElementById(anchorId);
        if (node.length > 0) {
            cy.animate({ fit: { eles: node, padding: 200 }, duration: 500 });
            node.addClass('highlighted-red');
            setTimeout(() => { node.removeClass('highlighted-red'); }, 1500);
        } else {
            alert('Anchor ID not found!');
        }
    });
    anchorIdInput.addEventListener('keyup', (e) => { if (e.key === 'Enter') viewAnchorBtn.click(); });

    dropAnchorBtn.addEventListener('click', () => {
        const inputText = anchorDropInput.value.trim();
        if (!inputText) return;

        const idsToDrop = inputText.replace(/,/g, ' ').split(/\s+/).filter(Boolean);

        // --- Validation Step ---
        const invalidIDs = idsToDrop.filter(id => !allAnchors.has(id));
        if (invalidIDs.length > 0) {
            alert(`The following anchor IDs were not found: ${invalidIDs.join(', ')}. No anchors were dropped.`);
            return;
        }

        const alreadyDropped = idsToDrop.filter(id => droppedAnchors.has(id));
        if (alreadyDropped.length > 0) {
            alert(`The following anchors have already been dropped: ${alreadyDropped.join(', ')}. No new anchors were dropped.`);
            return;
        }
        
        // --- Execution Step ---
        idsToDrop.forEach(id => {
            droppedAnchors.add(id);
            history.push({ type: 'anchor', id });
        });
        
        undoBtn.disabled = false;
        updateDroppedListsUI();
        updateGraph();
        anchorDropInput.value = '';
    });
    anchorDropInput.addEventListener('keyup', (e) => { if (e.key === 'Enter') dropAnchorBtn.click(); });


    dropReadBtn.addEventListener('click', () => {
        const inputText = readDropInput.value.trim();
        if (!inputText) return;

        const idsToDrop = inputText.replace(/,/g, ' ').split(/\s+/).filter(Boolean);
        
        // --- Validation Step ---
        const invalidIDs = idsToDrop.filter(id => !allReads.has(id));
        if (invalidIDs.length > 0) {
            alert(`The following read IDs were not found: ${invalidIDs.join(', ')}. No reads were dropped.`);
            return;
        }

        const alreadyDropped = idsToDrop.filter(id => droppedReads.has(id));
        if (alreadyDropped.length > 0) {
            alert(`The following reads have already been dropped: ${alreadyDropped.join(', ')}. No new reads were dropped.`);
            return;
        }
        
        // --- Execution Step ---
        idsToDrop.forEach(id => {
            droppedReads.add(id);
            history.push({ type: 'read', id });
        });
        
        undoBtn.disabled = false;
        updateDroppedListsUI();
        updateGraph();
        readDropInput.value = '';
    });
    readDropInput.addEventListener('keyup', (e) => { if (e.key === 'Enter') dropReadBtn.click(); });


    function handleRestore(e) {
        if (e.target.tagName === 'BUTTON') {
            const { type, id } = e.target.dataset;
            if (type === 'anchor') {
                droppedAnchors.delete(id);
            } else {
                droppedReads.delete(id);
            }
            // Restore is a direct action, does not affect undo history for drops
            updateDroppedListsUI();
            updateGraph();
        }
    }
    droppedAnchorsListDiv.addEventListener('click', handleRestore);
    droppedReadsListDiv.addEventListener('click', handleRestore);

    undoBtn.addEventListener('click', () => {
        if (history.length === 0) return;
        const lastAction = history.pop();
        if (lastAction.type === 'anchor') {
            droppedAnchors.delete(lastAction.id);
        } else {
            droppedReads.delete(lastAction.id);
        }
        undoBtn.disabled = history.length === 0;
        updateDroppedListsUI();
        updateGraph();
    });

    saveBtn.addEventListener('click', () => {
        const droppedSomething = droppedAnchors.size > 0 || droppedReads.size > 0;
        
        if (!droppedSomething && history.length === 0) {
             // A bit of a complex check: if user drops and then restores everything manually, nothing is left to save
            alert("No changes to save.");
            return;
        }

        if (!confirm("This will permanently modify the JSON files on the server. Are you sure?")) {
            return;
        }

        const dataToSave = {
            droppedAnchors: [...droppedAnchors],
            droppedReads: [...droppedReads]
        };

        fetch('/save', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(dataToSave)
        })
        .then(response => response.json())
        .then(result => {
            if (result.status === 'success') {
                alert('Changes saved successfully! The page will now reload.');
                location.reload();
            } else {
                alert(`Error saving changes: ${result.message}`);
            }
        })
        .catch(error => {
            console.error('Save error:', error);
            alert('An unexpected error occurred while trying to save.');
        });
    });

    // --- Helper Functions ---
    function updateStyles() {
        if (!cy) return;
        cy.style()
            .selector('node').style({ 'width': nodeSize, 'height': nodeSize })
            .selector('edge').style({ 'width': `mapData(weight, 1, 20, ${1 * edgeWidthMultiplier}, ${10 * edgeWidthMultiplier})`})
            .update();
    }

    function showNodeInfo(anchorId, reads) {
        let html = `<html><head><title>Anchor: ${anchorId}</title><style>body{font-family:sans-serif} table{width:100%; border-collapse:collapse} th,td{border:1px solid #ddd; padding:8px}</style></head><body>`;
        html += `<h3>Anchor: ${anchorId}</h3><h4>Reads:</h4>`;
        if (reads && reads.length > 0) {
            html += '<table><tr><th>Read ID</th><th>Strand</th><th>Start</th><th>End</th></tr>';
            reads.forEach(r => { html += `<tr><td>${r.read_id}</td><td>${r.read_strand}</td><td>${r.read_start}</td><td>${r.read_end}</td></tr>`; });
            html += '</table>';
        } else { html += '<p>No read information available.</p>'; }
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
            commonReads.forEach(r => { html += `<tr><td>${r.read_id}</td><td>${r.read_strand}</td><td>${r.read_start}</td><td>${r.read_end}</td></tr>`; });
            html += '</table>';
        } else { html += '<p>No common reads.</p>'; }
        html += '</body></html>';
        const newTab = window.open();
        newTab.document.write(html);
        newTab.document.close();
    }
    
    function resetHighlight() {
        if (!cy) return;
        cy.elements().removeClass('highlighted').removeClass('faded').removeClass('highlighted-red');
    }

    function highlightReadJourney(readId) {
        resetHighlight();
        
        const journey = originalReadInfo[readId]?.journey;

        if (!journey || journey.length === 0) {
            return; // Read has no defined journey
        }

        cy.elements().addClass('faded');

        // Highlight nodes in the journey
        const nodesToHighlight = cy.nodes().filter(node => journey.includes(node.id()));
        nodesToHighlight.removeClass('faded').addClass('highlighted');

        // Highlight edges between consecutive nodes in the journey
        for (let i = 0; i < journey.length - 1; i++) {
            const source = journey[i];
            const target = journey[i+1];

            // Check for the edge in both directions since the graph is undirected
            const edge1 = cy.edges(`[source = "${source}"][target = "${target}"]`);
            if (edge1.length > 0) {
                edge1.removeClass('faded').addClass('highlighted');
            }
            const edge2 = cy.edges(`[source = "${target}"][target = "${source}"]`);
            if (edge2.length > 0) {
                edge2.removeClass('faded').addClass('highlighted');
            }
        }
    }
});
