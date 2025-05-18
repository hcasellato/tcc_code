import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --------------------------
# 1. Data Loading with Validation
# --------------------------
def load_data(file_path):
    """Load data with validation checks"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Data file {file_path} not found!")
    
    try:
        df = pd.read_csv(file_path, sep=';', skiprows=2)
        required_columns = {'x', 'y', 'f(x,y)', 'Vx', 'Vy'}
        if not required_columns.issubset(df.columns):
            missing = required_columns - set(df.columns)
            raise ValueError(f"Missing columns: {missing}")
        return df
    except Exception as e:
        raise RuntimeError(f"Error loading {file_path}: {str(e)}")

# Load datasets with error handling
try:
    df_real = load_data("2D_finVolMet/true_2D_finVolMet_625.txt")
    df_test = load_data("2D_finVolMet/2D_finVolMet_625.txt")
except Exception as e:
    print(f"Fatal error: {str(e)}")
    exit(1)

# --------------------------
# 2. Heatmap Comparison
# --------------------------
def create_pivot(df):
    """Create pivot table with duplicate handling"""
    pivot = df.pivot_table(index='y', columns='x', values='f(x,y)', aggfunc='mean')
    if pivot.isna().any().any():
        print("Warning: NaN values detected in pivot table. Filling with zeros.")
        pivot = pivot.fillna(0)
    return pivot

# Create pivot tables
real_pivot = create_pivot(df_real)
test_pivot = create_pivot(df_test)

# Find common scale limits for color consistency
#vmin = min(real_pivot.min().min(), test_pivot.min().min())
#vmax = max(real_pivot.max().max(), test_pivot.max().max())

# Plot heatmaps
fig, ax = plt.subplots(1, 2, figsize=(18, 8), sharey=True)
cbar_ax = fig.add_axes([.91, .15, .02, .7])  # Shared colorbar axis

sns.heatmap(real_pivot, ax=ax[0], cbar_ax=cbar_ax,
            cmap='viridis')#, vmin=vmin, vmax=vmax)
ax[0].set_title('Real Data Heatmap')
ax[0].set_xlabel('X Coordinate')
ax[0].set_ylabel('Y Coordinate')

sns.heatmap(test_pivot, ax=ax[1], cbar_ax=cbar_ax, 
            cmap='viridis')#, vmin=vmin, vmax=vmax)
ax[1].set_title('Test Data Heatmap')
ax[1].set_xlabel('X Coordinate')

fig.suptitle("Scalar Field Comparison", y=0.95)
plt.savefig("comparacao_heatmaps.png", dpi=300, bbox_inches='tight')

# --------------------------
# 3. Vector Field Comparison (Modified Arrow Size)
# --------------------------
def plot_vector_field(ax, df, color, title, bol):
    """Plot vector field with optimized arrow scaling"""

    scale_value = 25 if bol == 1 else 100

    quiver = ax.quiver(
        df['x'], 
        df['y'], 
        df['Vx'], 
        df['Vy'], 
        color=color,
        angles='xy',
        scale_units='xy',
        scale = scale_value,
        width=0.002,  # Reduced shaft width
        headwidth=2.5,  # Reduced arrowhead size
        headlength=3, 
        headaxislength=2.5
    )
    ax.set_title(title)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.axis('equal')
    ax.set_xlim(0, 1)
    ax.set_ylim(0,1)
    return quiver

# Plot vector fields
fig, ax = plt.subplots(1, 2, figsize=(18, 8))
plot_vector_field(ax[0], df_real, 'blue', 'Real Vector Field',0)
plot_vector_field(ax[1], df_test, 'red', 'Test Vector Field',1)

fig.suptitle("Vector Field Comparison", y=0.95)
plt.tight_layout()
plt.savefig("comparacao_campos_vetoriais.png", dpi=300, bbox_inches='tight')