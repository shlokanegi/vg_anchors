import json
import os
import math

# --- Configuration ---
VISUALIZATION_DIR = 'visualization'
SNARLS_FILE = os.path.join(VISUALIZATION_DIR, 'snarl_to_anchors_dictionary.json')
READ_INFO_FILE = os.path.join(VISUALIZATION_DIR, 'read_info.json')
ANCHOR_READS_FILE = os.path.join(VISUALIZATION_DIR, 'anchors_read_info.json')
OUTPUT_FILE = os.path.join(VISUALIZATION_DIR, 'graph_layout.json')

# --- Layout Parameters ---
SNARL_SPACING = 250  # Horizontal space between snarls
INTRA_SNARL_SPACING = 30  # Space between anchors within a snarl
PULL_STRENGTH = 0.1  # Base strength for pulling distant nodes
ITERATIONS = 150 # Number of iterations for layout adjustment

# Initial layout waviness
WAVE_AMPLITUDE = 2000 # More vertical variation
WAVE_FREQUENCY = 0.08 # A bit wavier

# --- Helper Function ---
def get_snarl_numeric_id(snarl_id_str):
    """Extracts the first integer from a snarl ID string like 'snarl_id1-2-3'."""
    try:
        # Handle cases like 'snarl_id1' or '1-2-3'
        clean_id = snarl_id_str.replace('snarl_id', '')
        return int(clean_id.split('-')[0])
    except (ValueError, IndexError):
        return 0 # Fallback for unexpected formats

# --- Main Script ---
def generate_layout():
    print("Loading data files...")
    with open(SNARLS_FILE, 'r') as f:
        snarls_data = json.load(f)
    with open(READ_INFO_FILE, 'r') as f:
        read_info_data = json.load(f)
    with open(ANCHOR_READS_FILE, 'r') as f:
        anchor_reads_data = json.load(f)
    
    all_anchors = set(anchor_reads_data.keys())
    
    print("Calculating edge weights from read journeys...")
    edge_weights = {}
    for read_id, info in read_info_data.items():
        journey = info.get('journey', [])
        for i in range(len(journey) - 1):
            source, target = journey[i], journey[i+1]
            if source in all_anchors and target in all_anchors:
                key = tuple(sorted((source, target)))
                edge_weights[key] = edge_weights.get(key, 0) + 1

    # --- Initial Wavy Layout ---
    print("Performing initial wavy layout of snarls...")
    snarl_ids_sorted = sorted(snarls_data.keys(), key=get_snarl_numeric_id)
    node_positions = {}
    anchor_to_snarl_numeric = {}
    
    current_x = 0
    total_snarls = len(snarl_ids_sorted)
    for i, snarl_id in enumerate(snarl_ids_sorted):
        numeric_id = get_snarl_numeric_id(snarl_id)
        anchors_in_snarl = [a for a in snarls_data[snarl_id] if a in all_anchors]
        if not anchors_in_snarl:
            continue
        
        # Add a wave to the y-position to give it a curve
        progress = i / total_snarls if total_snarls > 0 else 0
        y_amplitude = WAVE_AMPLITUDE * math.sin(progress * math.pi) # Sine arch over the whole length
        snarl_y_center = y_amplitude * math.sin(i * WAVE_FREQUENCY)

        start_y_offset = - (len(anchors_in_snarl) - 1) * INTRA_SNARL_SPACING / 2
        for j, anchor_id in enumerate(anchors_in_snarl):
            node_positions[anchor_id] = {
                'x': current_x, 
                'y': snarl_y_center + start_y_offset + j * INTRA_SNARL_SPACING
            }
            anchor_to_snarl_numeric[anchor_id] = numeric_id
        
        current_x += SNARL_SPACING

    # --- Iterative Adjustment for Loopy Structure ---
    print(f"Running {ITERATIONS} iterations of layout adjustment...")
    for iteration in range(ITERATIONS):
        # The cooling factor reduces the strength of pulls over time, allowing layout to stabilize.
        cooling_factor = 1.0 - (iteration / ITERATIONS)

        for (source, target), weight in edge_weights.items():
            source_pos = node_positions.get(source)
            target_pos = node_positions.get(target)
            
            if not source_pos or not target_pos:
                continue

            source_snarl = anchor_to_snarl_numeric.get(source, 0)
            target_snarl = anchor_to_snarl_numeric.get(target, 0)

            snarl_dist = abs(source_snarl - target_snarl)

            # Only apply strong pull to distant snarls
            if snarl_dist > 1:
                # Dynamic pull strength based on edge weight and snarl distance
                dynamic_pull = PULL_STRENGTH * math.log1p(weight) * math.log1p(snarl_dist) * cooling_factor
                
                mid_x = (source_pos['x'] + target_pos['x']) / 2
                mid_y = (source_pos['y'] + target_pos['y']) / 2
                
                source_pos['x'] += (mid_x - source_pos['x']) * dynamic_pull
                source_pos['y'] += (mid_y - source_pos['y']) * dynamic_pull
                
                target_pos['x'] += (mid_x - target_pos['x']) * dynamic_pull
                target_pos['y'] += (mid_y - target_pos['y']) * dynamic_pull

    print("Preparing final output file...")
    output_graph = {
        'nodes': [],
        'edges': []
    }

    for anchor_id in all_anchors:
        if anchor_id in node_positions:
            output_graph['nodes'].append({
                'data': { 'id': anchor_id },
                'position': node_positions[anchor_id]
            })

    for (source, target), weight in edge_weights.items():
        output_graph['edges'].append({
            'data': {
                'source': source,
                'target': target,
                'weight': weight
            }
        })
    
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(output_graph, f, indent=4)

    print(f"Layout successfully generated and saved to {OUTPUT_FILE}")

if __name__ == '__main__':
    generate_layout() 
