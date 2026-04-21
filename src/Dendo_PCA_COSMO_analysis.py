import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import os
import sys
import warnings
import plotly.graph_objects as go

warnings.filterwarnings('ignore')

# ----------------------------------------------------------------------
# 1. Data loading and preprocessing
# ----------------------------------------------------------------------
def load_and_preprocess_data(file_path):
    """Load CSV, handle European decimals, extract features, names, ΔG."""
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

    # Try default decimal '.' first, fallback to ',' if needed
    try:
        df = pd.read_csv(file_path)
    except:
        try:
            df = pd.read_csv(file_path, decimal=',')
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            sys.exit(1)

    # Find compound column (Compounds or mol_name)
    comp_col = None
    for col in ['Compounds', 'mol_name']:
        if col in df.columns:
            comp_col = col
            break
    if comp_col is None:
        print("Error: No compound name column found ('Compounds' or 'mol_name').")
        sys.exit(1)

    # Find ΔG column
    dg_col = None
    for col in ['ΔG', 'dg']:
        if col in df.columns:
            dg_col = col
            break
    if dg_col is None:
        print("Error: No ΔG column found ('ΔG' or 'dg').")
        sys.exit(1)

    # Convert to numeric
    df[comp_col] = df[comp_col].astype(str)
    df[dg_col] = pd.to_numeric(df[dg_col], errors='coerce')

    # Descriptor columns: all numeric except comp_col and dg_col
    feature_cols = [c for c in df.columns if c not in [comp_col, dg_col]]
    for c in feature_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # Drop rows with missing critical values
    df_clean = df.dropna(subset=[comp_col, dg_col] + feature_cols).reset_index(drop=True)

    compound_names = df_clean[comp_col].values
    dg_values = df_clean[dg_col].values.astype(float)
    features = df_clean[feature_cols].values.astype(float)

    print(f"Loaded {len(compound_names)} compounds, {len(feature_cols)} descriptors.")
    return compound_names, dg_values, features, feature_cols

# ----------------------------------------------------------------------
# 2. Hierarchical clustering & BnB helpers
# ----------------------------------------------------------------------
def get_cluster_members(Z, cluster_id, n_samples):
    """Recursively return all leaf indices belonging to a cluster."""
    if cluster_id < n_samples:
        return {cluster_id}
    idx = int(cluster_id - n_samples)
    left, right = int(Z[idx, 0]), int(Z[idx, 1])
    return get_cluster_members(Z, left, n_samples) | get_cluster_members(Z, right, n_samples)

def find_medoid(indices, D):
    """Return index of the medoid (most central) in descriptor space."""
    if len(indices) == 1:
        return int(list(indices)[0])
    idx_list = list(indices)
    sub = D[np.ix_(idx_list, idx_list)]
    medoid_idx = idx_list[np.argmin(np.sum(sub, axis=1))]
    return int(medoid_idx)

def run_branch_and_bound(Z, n_samples, D, dg, max_depth=6, min_branch_size=2):
    """
    Branch‑and‑bound using dendrogram structure.
    Returns:
        results: list of terminal clusters with info
        rounds: list of dicts with round details (demasked compounds, kept branch members, etc.)
        stats: summary statistics
    """
    root = n_samples + len(Z) - 1
    demasked_so_far = set()
    rounds = []
    current_branches = [root]
    depth = 0
    terminal_results = []

    while current_branches:
        # Gather information for all current branches
        branch_info = []
        for cid in current_branches:
            members = get_cluster_members(Z, cid, n_samples)
            if not members:
                continue
            medoid_idx = find_medoid(members, D)
            branch_info.append((cid, members, medoid_idx))

        if not branch_info:
            break

        # Demask new medoids
        demasked_this_round = set()
        for _, _, medoid_idx in branch_info:
            if medoid_idx not in demasked_so_far:
                demasked_so_far.add(medoid_idx)
                demasked_this_round.add(medoid_idx)

        # Skip recording a round if no new demasking occurs (only possible when all medoids already known)
        if demasked_this_round:
            round_record = {
                'round_num': len(rounds) + 1,
                'demasked': demasked_this_round.copy(),
                'demasked_so_far': demasked_so_far.copy(),
            }
            rounds.append(round_record)

        # Compare branches by medoid ΔG
        with_dg = [(cid, members, medoid_idx, float(dg[medoid_idx])) for (cid, members, medoid_idx) in branch_info]
        with_dg.sort(key=lambda x: x[3])  # lower ΔG better

        # Separate terminal vs continue
        terminal, continue_b = [], []
        for cid, members, medoid_idx, medoid_dg in with_dg:
            if cid < n_samples or len(members) < min_branch_size or (max_depth is not None and depth >= max_depth):
                terminal.append((cid, members, medoid_idx, medoid_dg))
            else:
                continue_b.append((cid, members, medoid_idx, medoid_dg))
        terminal_results.extend(terminal)

        if not continue_b:
            if rounds:
                rounds[-1]['kept_branches'] = []
                rounds[-1]['filtered_branches'] = []
                rounds[-1]['kept_branch_members'] = {}
            break

        # Successive halving: keep best half
        n_keep = max(1, len(continue_b) // 2)
        kept = continue_b[:n_keep]
        filtered = continue_b[n_keep:]

        if rounds:
            rounds[-1]['kept_branches'] = [b[0] for b in kept]
            rounds[-1]['filtered_branches'] = [b[0] for b in filtered]
            rounds[-1]['kept_branch_members'] = {b[0]: b[1] for b in kept}

        # Expand children for next round
        next_branches = []
        for cid, _, _, _ in kept:
            if cid < n_samples:
                next_branches.append(cid)
            else:
                idx = cid - n_samples
                next_branches.append(int(Z[idx, 0]))
                next_branches.append(int(Z[idx, 1]))
        current_branches = next_branches
        depth += 1

    # Summary statistics
    stats = {
        'total_rounds': len(rounds),
        'total_compounds_demasked': len(demasked_so_far),
        'global_best_dg': np.min(dg),
        'best_dg_found': None
    }
    if terminal_results:
        all_terminal_members = set()
        for _, members, _, _ in terminal_results:
            all_terminal_members.update(members)
        best_dg_found = np.min(dg[list(all_terminal_members)])
        stats['best_dg_found'] = best_dg_found
    else:
        stats['best_dg_found'] = np.inf

    return terminal_results, rounds, stats

# ----------------------------------------------------------------------
# 3. Dendrogram plotting (Matplotlib + Plotly) with demasked highlighting
# ----------------------------------------------------------------------
def plot_dendrogram_with_dg(Z, compounds, dg_values, demasked_indices, output_base):
    """
    Create a vertical dendrogram with leaves coloured by ΔG.
    Leaves that were never demasked are shown in grey.
    """
    # Matplotlib version for PNG/SVG
    fig, ax = plt.subplots(figsize=(16, 10))
    dendro = dendrogram(
        Z,
        labels=compounds.tolist(),
        orientation='top',
        ax=ax,
        leaf_font_size=7,
        color_threshold=None,
    )
    leaf_order = dendro['leaves']
    n_leaves = len(leaf_order)
    leaf_x = [5 + i*10 for i in range(n_leaves)]
    leaf_y = [0.0] * n_leaves

    dg_ordered = dg_values[leaf_order]
    demasked_set = set(demasked_indices)
    is_demasked = np.array([i in demasked_set for i in leaf_order])

    cmap = plt.get_cmap('RdYlGn_r')
    demasked_dg = dg_ordered[is_demasked]
    if len(demasked_dg) > 0:
        norm = mpl.colors.Normalize(vmin=np.min(demasked_dg), vmax=np.max(demasked_dg))
    else:
        norm = mpl.colors.Normalize(vmin=np.min(dg_ordered), vmax=np.max(dg_ordered))

    colors = []
    for i, dg_val in enumerate(dg_ordered):
        if is_demasked[i]:
            colors.append(cmap(norm(dg_val)))
        else:
            colors.append('white')

    scatter = ax.scatter(leaf_x, leaf_y, c=colors, s=150,
                         edgecolors='black', linewidths=0.8, zorder=10)

    if len(demasked_dg) > 0:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_label('ΔG (kcal/mol)', rotation=270, labelpad=15)

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin - 0.05*(ymax-ymin), ymax)
    ax.set_title('Hierarchical Clustering (Ward, Euclidean)\nLeaves coloured by ΔG (demasked only)',
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('Compound')
    ax.set_ylabel('Distance')
    plt.tight_layout()
    fig.savefig(f"{output_base}.png", dpi=150, bbox_inches='tight')
    fig.savefig(f"{output_base}.svg", format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"Saved dendrogram: {output_base}.png / .svg")
    
# ----------------------------------------------------------------------
# 4. Main
# ----------------------------------------------------------------------
def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <csv_file> <max_depth> <min_branch_size>")
        print("Example: python dendro_analysis.py data.csv 6 2")
        sys.exit(1)

    input_file = sys.argv[1]
    try:
        max_depth = int(sys.argv[2])
        min_branch_size = int(sys.argv[3])
    except ValueError:
        print("Error: max_depth and min_branch_size must be integers.")
        sys.exit(1)

    print(f"Input file: {input_file}")
    print(f"BnB parameters: max_depth={max_depth}, min_branch_size={min_branch_size}")

    # Load data
    compounds, dg_values, features, feature_cols = load_and_preprocess_data(input_file)

    # Standardise
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(features)

    # Hierarchical clustering (Ward, Euclidean)
    print("Performing hierarchical clustering (Ward, Euclidean)...")
    Z = linkage(X_scaled, method='ward', metric='euclidean')
    D = squareform(pdist(X_scaled, metric='euclidean'))

    # Run BnB
    print("Running Branch‑and‑Bound search...")
    terminal_results, rounds, stats = run_branch_and_bound(
        Z, len(compounds), D, dg_values,
        max_depth=max_depth, min_branch_size=min_branch_size
    )

    global_best_dg = stats['global_best_dg']
    global_best_idx = np.argmin(dg_values)
    best_dg_found = stats['best_dg_found']

    print(f"BnB completed: {stats['total_rounds']} rounds, {stats['total_compounds_demasked']} compounds demasked.")
    print(f"Best ΔG found: {best_dg_found:.3f} (global best: {global_best_dg:.3f})")

    # Collect all demasked indices for plotting
    demasked_all = set()
    for rd in rounds:
        demasked_all.update(rd['demasked'])

    # Generate dendrogram plots
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    plot_base = f"{base_name}_dendrogram"
    plot_dendrogram_with_dg(Z, compounds, dg_values, demasked_all, plot_base)

    # Write rounds assignment file
    rounds_file = f"{base_name}_rounds_assignment.txt"
    with open(rounds_file, 'w') as f:
        f.write("Branch‑and‑Bound Rounds Assignment\n")
        f.write("="*60 + "\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"Clustering: Ward, Euclidean\n")
        f.write(f"BnB parameters: max_depth={max_depth}, min_branch_size={min_branch_size}\n")
        f.write(f"Total compounds: {len(compounds)}\n")
        f.write(f"Global best ΔG: {global_best_dg:.3f} (Compound: {compounds[global_best_idx]})\n")
        f.write(f"Best ΔG found by BnB: {best_dg_found:.3f}\n")
        f.write("="*60 + "\n\n")

        for rd in rounds:
            round_num = rd['round_num']
            demasked_names = [compounds[i] for i in rd['demasked']]
            f.write(f"Round {round_num}: demasked {len(demasked_names)} compound(s): {', '.join(demasked_names)}\n")
            if 'kept_branches' in rd:
                kept_ids = rd['kept_branches']
                filtered_ids = rd.get('filtered_branches', [])
                f.write(f"  Kept branches: {kept_ids}\n")
                f.write(f"  Filtered branches: {filtered_ids}\n")
                for cid in kept_ids:
                    members = rd['kept_branch_members'].get(cid, set())
                    member_names = [compounds[i] for i in members]
                    f.write(f"    Branch {cid} members: {', '.join(member_names)}\n")
            f.write("\n")

        f.write("Terminal clusters after BnB:\n")
        for cid, members, medoid_idx, medoid_dg in terminal_results:
            member_names = [compounds[i] for i in members]
            f.write(f"  Cluster {cid}: medoid {compounds[medoid_idx]} (ΔG={medoid_dg:.3f}), members: {', '.join(member_names)}\n")

    print(f"Rounds assignment saved to: {rounds_file}")

if __name__ == "__main__":
    main()
