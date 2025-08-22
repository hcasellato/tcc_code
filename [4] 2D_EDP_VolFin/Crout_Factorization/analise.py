import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --------------------------
# 1. Data Loading with Validation
# --------------------------
def load_data(file_path):
    """Load data with validation checks
    
    Args:
        file_path (str): Path to the data file
        
    Returns:
        pd.DataFrame: Loaded dataframe
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns are missing
        RuntimeError: For other loading errors
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Data file {file_path} not found!")
    
    try:
        df = pd.read_csv(file_path, sep=';', skiprows=2)
        required_columns = {'x', 'y', 'f(x,y)', 'exact_solution', 'difference', 'Vx', 'Vy', 'VxE', 'VyE'}
        if not required_columns.issubset(df.columns):
            missing = required_columns - set(df.columns)
            raise ValueError(f"Missing columns: {missing}")
        return df
    except Exception as e:
        raise RuntimeError(f"Error loading {file_path}: {str(e)}")

# --------------------------
# 2. Heatmap Comparison
# --------------------------
def plot_heatmaps(df, output_path="comparacao_heatmaps.png"):
    """Plot comparison heatmaps for scalar fields
    
    Args:
        df (pd.DataFrame): Input data
        output_path (str): Path to save the figure
    """
    try:
        # Create pivot tables
        pivot_test = df.pivot_table(index='y', columns='x', values='f(x,y)', aggfunc='mean')
        pivot_real = df.pivot_table(index='y', columns='x', values='exact_solution', aggfunc='mean')
        pivot_diff = df.pivot_table(index='y', columns='x', values='difference', aggfunc='mean')

        # Create figure with 3 subplots
        fig, ax = plt.subplots(1, 3, figsize=(24, 8), sharey=True)
        fig.suptitle("Scalar Field Comparison", y=0.95)
        
        # Plot heatmaps
        sns.heatmap(pivot_real, ax=ax[0], cmap='viridis', vmin=-1, vmax=1)
        ax[0].set_title('Exact Solution')
        ax[0].set_xlabel('X Coordinate')
        ax[0].set_ylabel('Y Coordinate')

        sns.heatmap(pivot_test, ax=ax[1], cmap='viridis', vmin=-1, vmax=1)
        ax[1].set_title('Numerical Solution')
        ax[1].set_xlabel('X Coordinate')

        sns.heatmap(pivot_diff, ax=ax[2], cmap='coolwarm', vmin=-1, vmax=1)
        ax[2].set_title('Difference')
        ax[2].set_xlabel('X Coordinate')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error plotting heatmaps: {str(e)}")
        raise

# --------------------------
# 3. Vector Field Comparison
# --------------------------
def plot_vector_fields(df, output_path="comparacao_campos_vetoriais.png"):
    """Plot comparison of vector fields
    
    Args:
        df (pd.DataFrame): Input data containing vector fields
        output_path (str): Path to save the figure
    """
    try:
        fig, ax = plt.subplots(1, 2, figsize=(18, 8))
        plt.style.use('seaborn-v0_8-pastel')
        
        # Plot exact vector field
        ax[0].quiver(
            df['x'], df['y'], 
            df['VxE'], df['VyE'],
            color='blue', angles='xy', scale_units='xy',
            scale=50, width=0.002, headwidth=2.5
        )
        ax[0].set_title('Exact Vector Field')
        ax[0].set_xlabel('X Coordinate')
        ax[0].set_ylabel('Y Coordinate')
        
        # Plot numerical vector field
        ax[1].quiver(
            df['x'], df['y'], 
            df['Vx'], df['Vy'],
            color='red', angles='xy', scale_units='xy',
            scale=50, width=0.002, headwidth=2.5
        )
        ax[1].set_title('Numerical Vector Field')
        ax[1].set_xlabel('X Coordinate')
        
        # Set common properties
        for a in ax:
            a.axis('equal')
            a.set_xlim(df['x'].min(), df['x'].max())
            a.set_ylim(df['y'].min(), df['y'].max())
        
        fig.suptitle("Vector Field Comparison", y=0.95)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error plotting vector fields: {str(e)}")
        raise

# --------------------------
# Main Execution
# --------------------------
if __name__ == "__main__":
    try:
        # Load data
        df = load_data("2D_finVolMet/2D_finVolMet_625.txt")
        
        # Generate visualizations
        plot_heatmaps(df)
        plot_vector_fields(df)
        
    except Exception as e:
        print(f"Execution failed: {str(e)}")
        exit(1)