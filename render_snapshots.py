import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio

def read_ply(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    start_index = 0
    for i, line in enumerate(lines):
        if line.strip() == 'end_header':
            start_index = i + 1
            break
    
    pts = []
    norms = []
    ids = []
    for i in range(start_index, len(lines)):
        parts = lines[i].split()
        if len(parts) >= 7:
            pts.append([float(parts[0]), float(parts[1]), float(parts[2])])
            norms.append([float(parts[3]), float(parts[4]), float(parts[5])])
            ids.append(int(parts[6]))
    return np.array(pts), np.array(norms), np.array(ids)

def render_snapshot(ply_path, output_path, gen, bounds=None):
    pts, norms, ids = read_ply(ply_path)
    if len(pts) == 0:
        return
    
    fig = plt.figure(figsize=(10, 10), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')
    
    # --- 3D Shading Simulation ---
    # Virtual light source from top-right
    light_dir = np.array([1.0, 1.0, 2.0])
    light_dir /= np.linalg.norm(light_dir)
    
    # Calculate shading intensity (Lambertian)
    shading = np.sum(norms * light_dir, axis=1)
    shading = np.clip(shading, 0.2, 1.0) # Ambient + Diffuse
    
    # --- Categorical Coloring ---
    cmap = plt.cm.get_cmap('tab20')
    unique_ids = np.unique(ids)
    
    for sid in unique_ids:
        mask = ids == sid
        shard_pts = pts[mask]
        shard_shading = shading[mask]
        
        # Base color for this shard
        base_color = np.array(cmap(sid % 20)[:3])
        # Apply shading to base color
        colors = base_color * shard_shading[:, np.newaxis]
        
        # Render with larger points for "solid" look
        ax.scatter(shard_pts[:, 0], shard_pts[:, 1], shard_pts[:, 2], 
                   s=15, c=colors, alpha=1.0, edgecolors='none')

    ax.set_title(f'Pottery Reconstruction - Gen {gen}', fontsize=22, pad=-30, color='#3E2723', fontweight='bold')
    
    # --- Auto-Orientation & Scaling ---
    if bounds is not None:
        ax.set_xlim(bounds[0], bounds[1])
        ax.set_ylim(bounds[2], bounds[3])
        ax.set_zlim(bounds[4], bounds[5])
    else:
        # Centering on current points
        center = pts.mean(axis=0)
        max_range = (pts.max(axis=0) - pts.min(axis=0)).max() / 2.0
        ax.set_xlim(center[0] - max_range, center[0] + max_range)
        ax.set_ylim(center[1] - max_range, center[1] + max_range)
        ax.set_zlim(center[2] - max_range, center[2] + max_range)

    # Clean UI
    ax.set_axis_off()
    ax.grid(False)
    
    # Best angle for viewing pots (slightly from above)
    ax.view_init(elev=20, azim=-45)
    
    plt.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close()

def main():
    result_dir = 'result_paper'
    snapshot_files = sorted(glob.glob(os.path.join(result_dir, 'ga_snapshot_gen_*.ply')))
    
    if not snapshot_files:
        print("No snapshots found. Run C++ first.")
        return

    # Calculate global bounds for stability in GIF
    all_pts, _, _ = read_ply(snapshot_files[-1]) # Use last one for max extent
    center = all_pts.mean(axis=0)
    max_range = (all_pts.max(axis=0) - all_pts.min(axis=0)).max() / 2.0
    bounds = [center[0]-max_range, center[0]+max_range, 
              center[1]-max_range, center[1]+max_range, 
              center[2]-max_range, center[2]+max_range]

    frames = []
    print(f"Found {len(snapshot_files)} snapshots. Rendering with shading...")

    for ply_file in snapshot_files:
        gen_str = ply_file.split('_')[-1].replace('.ply', '')
        gen = int(gen_str)
        img_path = ply_file.replace('.ply', '.png')
        
        print(f" Rendering Gen {gen}...")
        render_snapshot(ply_file, img_path, gen, bounds)
        frames.append(imageio.imread(img_path))

    # Create GIF
    gif_path = os.path.join(result_dir, 'assembly_evolution.gif')
    imageio.mimsave(gif_path, frames, fps=3)
    print(f"Enhanced GIF saved to {gif_path}")

if __name__ == "__main__":
    main()
