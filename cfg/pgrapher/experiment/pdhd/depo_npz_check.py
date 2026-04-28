import numpy as np
data = np.load('depositions.npz')

# Check exact shapes and types
print("depo_data_0:")
print(f"  Shape: {data['depo_data_0'].shape}")
print(f"  Columns: {data['depo_data_0'].shape[1] if len(data['depo_data_0'].shape) > 1 else 'NOT 2D'}")

print("depo_info_0:")  
print(f"  Shape: {data['depo_info_0'].shape}")
print(f"  Columns: {data['depo_info_0'].shape[1] if len(data['depo_info_0'].shape) > 1 else 'NOT 2D'}")

# Check if rows match
print(f"Rows match: {data['depo_data_0'].shape[0] == data['depo_info_0'].shape[0]}")