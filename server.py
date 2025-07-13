import os
import json
from flask import Flask, request, jsonify, send_from_directory

app = Flask(__name__)

# The directory where your visualization files (html, css, js, json) are located
VISUALIZATION_DIR = 'visualization'

# --- API Endpoint to Save Changes ---
@app.route('/save', methods=['POST'])
def save_changes():
    try:
        data = request.get_json()
        dropped_anchors = data.get('droppedAnchors', [])
        dropped_reads = data.get('droppedReads', [])

        if not dropped_anchors and not dropped_reads:
            return jsonify({"status": "nothing to save"}), 200

        # --- Update snarl_to_anchors_dictionary.json ---
        snarl_path = os.path.join(VISUALIZATION_DIR, 'snarl_to_anchors_dictionary.json')
        with open(snarl_path, 'r+') as f:
            snarl_data = json.load(f)
            for snarl_id in snarl_data:
                snarl_data[snarl_id] = [anchor for anchor in snarl_data[snarl_id] if anchor not in dropped_anchors]
            f.seek(0)
            f.truncate()
            json.dump(snarl_data, f, indent=4)

        # --- Update anchors_read_info.json ---
        anchor_reads_path = os.path.join(VISUALIZATION_DIR, 'anchors_read_info.json')
        with open(anchor_reads_path, 'r+') as f:
            anchor_reads_data = json.load(f)
            # Remove dropped anchors
            anchor_reads_data = {anchor: anchor_reads_data[anchor] for anchor in anchor_reads_data if anchor not in dropped_anchors}
            # Remove dropped reads from remaining anchors
            for anchor in anchor_reads_data:
                anchor_reads_data[anchor] = [read for read in anchor_reads_data[anchor] if read['read_id'] not in dropped_reads]
            f.seek(0)
            f.truncate()
            json.dump(anchor_reads_data, f, indent=4)

        # # --- Update read_info.json ---
        # read_info_path = os.path.join(VISUALIZATION_DIR, 'read_info.json')
        # with open(read_info_path, 'r+') as f:
        #     read_info_data = json.load(f)
        #     read_info_data = {read: read_info_data[read] for read in read_info_data if read not in dropped_reads}
        #     f.seek(0)
        #     f.truncate()
        #     json.dump(read_info_data, f, indent=4)
            
        return jsonify({"status": "success"}), 200

    except Exception as e:
        print(f"Error saving changes: {e}")
        return jsonify({"status": "error", "message": str(e)}), 500


# --- Static File Serving ---
@app.route('/')
def root():
    return send_from_directory(VISUALIZATION_DIR, 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    # This is a bit of a catch-all. It serves files from the visualization directory.
    return send_from_directory(VISUALIZATION_DIR, path)

if __name__ == '__main__':
    app.run(debug=True, port=8000) 
