from assembler.anchor import Anchor
from collections import Counter
import assembler.config as settings


def find_adjacent_snarl_pairs(snarl_list, snarl_to_anchors_dictionary, read_to_snarl_dictionary):
    """
    This function uses the read journeys to find adjacent snarl pairs.
    """
    adjacent_snarl_pairs = []

    for current_snarl_id in snarl_list:
        current_snarl_anchors = snarl_to_anchors_dictionary[current_snarl_id]
        # print(f"current snarl ID: {current_snarl_id} has {len(current_snarl_anchors)} anchors", flush=True)

        potential_prev_snarl_ids = Counter()
        potential_succ_snarl_ids = Counter()        
        
        for anchor in current_snarl_anchors:
            # print(f" anchor: {anchor!r} has {len(anchor.bp_matched_reads)} reads", flush=True)
            for read_id, strand, *_ in anchor.bp_matched_reads:
                anchor_rank = anchor.read_ranks[read_id]
                read_journey_snarls = read_to_snarl_dictionary[read_id]
                # print(f"  read ID: {read_id}, strand: {strand}, anchor rank: {anchor_rank}", flush=True)
                # print(f"  read journey snarls: {read_journey_snarls}", flush=True)
                
                # snarl before the current anchor
                prev_snarl = (
                    read_journey_snarls[anchor_rank - 1] if anchor_rank > 0 else "NA"
                )

                # snarl after the current anchor
                succ_snarl = (
                    read_journey_snarls[anchor_rank + 1] if anchor_rank + 1 < len(read_journey_snarls) else "NA"
                )

                if strand == 0:
                    potential_prev_snarl_ids[prev_snarl] += 1
                    potential_succ_snarl_ids[succ_snarl] += 1
                    # print(f"  prev snarl: {prev_snarl}", flush=True)
                    # print(f"  succ snarl: {succ_snarl}", flush=True)
                else:
                    potential_prev_snarl_ids[succ_snarl] += 1
                    potential_succ_snarl_ids[prev_snarl] += 1
                    # print(f"  prev snarl: {succ_snarl}", flush=True)
                    # print(f"  succ snarl: {prev_snarl}", flush=True)

        if not potential_prev_snarl_ids and not potential_succ_snarl_ids:
            continue

        # print(f"   potential prev snarls: {potential_prev_snarl_ids}", flush=True)
        # print(f"   potential succ snarls: {potential_succ_snarl_ids}", flush=True)

        # select the most frequent prev and succ snarls
        most_frequent_prev_snarl_id = (
            potential_prev_snarl_ids.most_common(1)[0][0] if potential_prev_snarl_ids else "NA"
        )
        most_frequent_succ_snarl_id = (
            potential_succ_snarl_ids.most_common(1)[0][0] if potential_succ_snarl_ids else "NA"
        )
    
        if most_frequent_prev_snarl_id != "NA":
            adjacent_snarl_pairs.append((most_frequent_prev_snarl_id, current_snarl_id))
        if most_frequent_succ_snarl_id != "NA":
            adjacent_snarl_pairs.append((current_snarl_id, most_frequent_succ_snarl_id))
        
        # print(f"   adjacent snarl pairs: {adjacent_snarl_pairs}", flush=True)

    return adjacent_snarl_pairs
