import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA
import plotly.graph_objects as go
import plotly.io as pio
from collections import defaultdict
import warnings
import sys
import os
warnings.filterwarnings('ignore')

def load_and_preprocess_data(file_path):
    """Load and preprocess the QTAIM data"""
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    
    try:
        # Read the CSV file with European decimal format (comma as decimal separator)
        df = pd.read_csv(file_path, decimal=',')
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)
    
    # Check if 'Compounds' column exists
    if 'Compounds' not in df.columns:
        print("Error: 'Compounds' column not found in the CSV file.")
        sys.exit(1)
    
    # Extract feature columns (all columns except 'Compounds')
    feature_columns = [col for col in df.columns if col != 'Compounds']
    
    # Separate compound names and features
    compound_names = df['Compounds'].values
    features = df[feature_columns].values
    
    return compound_names, features, feature_columns

def perform_scaling_and_pca(features, n_components=2, random_state=42):
    """Perform MinMax scaling and PCA transformation"""
    # MinMax scaling
    scaler = MinMaxScaler()
    scaled_features = scaler.fit_transform(features)
    
    # PCA transformation
    pca = PCA(n_components=n_components, random_state=random_state)
    pca_results = pca.fit_transform(scaled_features)
    
    # Calculate variance explained by each component
    explained_variance = pca.explained_variance_ratio_ * 100
    
    return pca_results, explained_variance

def assign_to_grid(pca_results, n_rows, m_cols):
    """Assign each point to its grid cell based on user-defined dimensions"""
    # Normalize PCA coordinates to [0, 1] range
    pca_normalized = (pca_results - pca_results.min(axis=0)) / (pca_results.max(axis=0) - pca_results.min(axis=0))
    
    # Calculate grid coordinates
    grid_x = (pca_normalized[:, 0] * m_cols).astype(int)  # X corresponds to columns
    grid_y = (pca_normalized[:, 1] * n_rows).astype(int)  # Y corresponds to rows
    
    # Ensure coordinates are within bounds
    grid_x = np.clip(grid_x, 0, m_cols - 1)
    grid_y = np.clip(grid_y, 0, n_rows - 1)
    
    # Create grid assignment dictionary
    grid_assignments = {}
    cell_contents = defaultdict(list)
    
    for i, (x, y) in enumerate(zip(grid_x, grid_y)):
        cell_id = (x, y)  # (column, row)
        grid_assignments[i] = cell_id
        cell_contents[cell_id].append(i)
    
    return grid_assignments, cell_contents

def generate_rounds(compound_names, cell_contents):
    """Generate rounds with one entry from each populated grid cell"""
    # Create a copy of cell contents to track remaining entries
    remaining_entries = {cell: entries.copy() for cell, entries in cell_contents.items()}
    
    rounds = []
    round_num = 1
    
    while any(remaining_entries.values()):
        current_round = []
        
        # For each cell that still has entries, take one
        for cell in list(remaining_entries.keys()):
            if remaining_entries[cell]:
                entry_idx = remaining_entries[cell].pop(0)
                current_round.append((entry_idx, compound_names[entry_idx]))
                
                # Remove empty cells
                if not remaining_entries[cell]:
                    del remaining_entries[cell]
        
        # Sort by index to maintain order
        current_round.sort(key=lambda x: x[0])
        rounds.append((round_num, [name for idx, name in current_round]))
        round_num += 1
    
    return rounds

def plot_pca_results_plotly(pca_results, compound_names, rounds, cell_contents, n_rows, m_cols, output_file='pca_visualization.html'):
    """Create an interactive Plotly visualization of PCA results with grid and colored rounds"""
    # Create color palette for rounds - using high contrast colors
    import plotly.express as px
    colors = px.colors.qualitative.Bold
    
    # Create a mapping from compound index to round number and color
    compound_to_round = {}
    compound_to_color = {}
    
    for round_num, round_entries in rounds:
        round_color = colors[(round_num - 1) % len(colors)]
        for name in round_entries:
            idx = np.where(compound_names == name)[0][0]
            compound_to_round[idx] = round_num
            compound_to_color[idx] = round_color
    
    # Create scatter plot traces
    traces = []
    
    # Group compounds by round for better visualization
    for round_num in sorted(set(compound_to_round.values())):
        round_indices = [idx for idx, rnd in compound_to_round.items() if rnd == round_num]
        round_compounds = [compound_names[idx] for idx in round_indices]
        round_color = compound_to_color[round_indices[0]]
        
        trace = go.Scatter(
            x=pca_results[round_indices, 0],
            y=pca_results[round_indices, 1],
            mode='markers+text',
            name=f'Round {round_num}',
            text=round_compounds,
            textposition="top center",
            marker=dict(
                size=20,  # Increased marker size
                color=round_color,
                line=dict(width=3, color='black')  # Thicker black border
            ),
            textfont=dict(size=16, family='Arial', weight='bold', color='black'),  # Increased font size, black text
            hovertemplate=(
                '<b>%{text}</b><br>' +
                'PC1: %{x:.3f}<br>' +
                'PC2: %{y:.3f}<br>' +
                f'Round: {round_num}<br>' +
                '<extra></extra>'
            )
        )
        traces.append(trace)
    
    # Create figure
    fig = go.Figure(data=traces)
    
    # Add grid lines - thicker and darker for better visibility
    x_min, x_max = pca_results[:, 0].min(), pca_results[:, 0].max()
    y_min, y_max = pca_results[:, 1].min(), pca_results[:, 1].max()
    
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    # Add vertical grid lines
    for i in range(1, m_cols):
        x_val = x_min + (i / m_cols) * x_range
        fig.add_shape(
            type="line",
            x0=x_val, y0=y_min,
            x1=x_val, y1=y_max,
            line=dict(color="gray", width=2, dash="dash"),  # Thicker grid lines
            layer="below"
        )
    
    # Add horizontal grid lines
    for i in range(1, n_rows):
        y_val = y_min + (i / n_rows) * y_range
        fig.add_shape(
            type="line",
            x0=x_min, y0=y_val,
            x1=x_max, y1=y_val,
            line=dict(color="gray", width=2, dash="dash"),  # Thicker grid lines
            layer="below"
        )
    
    
    # Update layout with high visibility settings
    fig.update_layout(
        title=dict(
            text=f'PCA Visualization with {n_rows}x{m_cols} Grid (Colored by Round)',
            font=dict(size=28, family='Arial', weight='bold', color='black'),  # Larger title
            x=0.5,
            xanchor='right'
        ),
        xaxis=dict(
            title=dict(text='PC1', font=dict(size=24, weight='bold', color='black')),  # Larger axis label
            gridcolor='lightgray',
            showgrid=False,
            tickfont=dict(size=20, color='black'),  # Larger tick labels
            linecolor='black',  # Axis line color
            linewidth=3,  # Thicker axis line
            mirror=False,
            showline=True,
            zeroline=False  # Don't show default zero line since we added our own
        ),
        yaxis=dict(
            title=dict(text='PC2', font=dict(size=24, weight='bold', color='black')),  # Larger axis label
            gridcolor='lightgray',
            showgrid=False,
            tickfont=dict(size=20, color='black'),  # Larger tick labels
            linecolor='black',  # Axis line color
            linewidth=3,  # Thicker axis line
            mirror=False,
            showline=True,
            zeroline=False  # Don't show default zero line since we added our own
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        showlegend=True,
        legend=dict(
            title=dict(text='Rounds', font=dict(size=22, weight='bold', color='black')),  # Larger legend title
            font=dict(size=20, color='black'),  # Larger legend text
            x=1.05,
            y=1,
            xanchor='left',
            yanchor='top',
            #bordercolor='black',
            #borderwidth=2,
            #bgcolor='rgba(255, 255, 255, 0.9)'
        ),
        width=1400,  # Increased width
        height=1000,  # Increased height
        hoverlabel=dict(
            bgcolor="white",
            font_size=18,  # Larger hover text
            font_family="Arial",
            font_color="black"
        ),
        margin=dict(l=100, r=200, t=100, b=100)  # Adjust margins for larger elements
    )
    
    # Save as HTML (interactive) and PNG (static)
    html_file = output_file
    png_file = output_file.replace('.html', '.png')
    
    fig.write_html(html_file)
    print(f"Interactive plot saved as: {html_file}")
    
    # Also save as PNG for static version
    fig.write_image(png_file, width=1400, height=1000, scale=2)  # Higher resolution
    print(f"Static plot saved as: {png_file}")
    
    return fig

def main():
    # Parse command line arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <csv_file> <n_rows> <m_cols>")
        print("Example: python Bader_PCA_analysis.py QTAIM_1-24_noTS_ver3.csv 5 5")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    try:
        n_rows = int(sys.argv[2])
        m_cols = int(sys.argv[3])
    except ValueError:
        print("Error: n_rows and m_cols must be integers.")
        sys.exit(1)
    
    # Validate grid dimensions
    if n_rows <= 0 or m_cols <= 0:
        print("Error: Grid dimensions must be positive integers.")
        sys.exit(1)
    
    print(f"Input file: {input_file}")
    print(f"Grid dimensions: {n_rows} rows x {m_cols} columns")
    
    # Generate output filenames based on input filename
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}_rounds_assignment.txt"
    plot_file = f"{base_name}_pca_visualization.html"
    
    print("Loading and preprocessing data...")
    compound_names, features, feature_columns = load_and_preprocess_data(input_file)
    
    print("Performing MinMax scaling and PCA transformation...")
    pca_results, explained_variance = perform_scaling_and_pca(features)
    
    print(f"PCA explained variance: PC1: {explained_variance[0]:.2f}%, PC2: {explained_variance[1]:.2f}%")
    print(f"Using grid size: {n_rows} rows x {m_cols} columns")
    
    print("Assigning compounds to grid cells...")
    grid_assignments, cell_contents = assign_to_grid(pca_results, n_rows, m_cols)
    
    # Print grid statistics
    cell_counts = [len(entries) for entries in cell_contents.values()]
    if cell_counts:
        mean_entries = np.mean(cell_counts)
        print(f"Grid statistics: {len(cell_contents)} populated cells out of {n_rows * m_cols} total")
        print(f"Min entries per cell: {min(cell_counts)}")
        print(f"Max entries per cell: {max(cell_counts)}")
        print(f"Mean entries per cell: {mean_entries:.2f}")
    else:
        print("No cells populated!")
    
    print("Generating rounds...")
    rounds = generate_rounds(compound_names, cell_contents)
    
    # Create visualization with Plotly
    print(f"Creating visualization...")
    plot_pca_results_plotly(pca_results, compound_names, rounds, cell_contents, n_rows, m_cols, plot_file)
    
    # Write results to file
    print(f"Writing results to {output_file}...")
    with open(output_file, 'w') as f:
        # Write header with variance information
        f.write(f"QTAIM Data Rounds Assignment\n")
        f.write("=" * 50 + "\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"PCA Explained Variance:\n")
        f.write(f"  PC1: {explained_variance[0]:.2f}%\n")
        f.write(f"  PC2: {explained_variance[1]:.2f}%\n")
        f.write(f"Grid size: {n_rows} rows x {m_cols} columns\n")
        f.write(f"Total compounds: {len(compound_names)}\n")
        f.write("=" * 50 + "\n\n")
        
        # Write grid information
        f.write("Grid Cell Contents (format: (column, row)):\n")
        for cell, entries in sorted(cell_contents.items()):
            entry_names = [compound_names[idx] for idx in entries]
            f.write(f"  Cell {cell}: {', '.join(entry_names)}\n")
        f.write("\n")
        
        # Write rounds
        f.write("Rounds Assignment:\n")
        for round_num, round_entries in rounds:
            f.write(f"Round {round_num}: {', '.join(round_entries)}\n")
    
    print("Process completed successfully!")
    print(f"Results saved to: {output_file}")
    print(f"Interactive plot saved as: {plot_file}")
    print(f"Static PNG plot saved as: {plot_file.replace('.html', '.png')}")
    
    # Print summary to console
    print("\nSummary:")
    print(f"Total rounds generated: {len(rounds)}")
    for round_num, round_entries in rounds:
        print(f"Round {round_num}: {len(round_entries)} compounds - {', '.join(round_entries)}")

if __name__ == "__main__":
    main()
