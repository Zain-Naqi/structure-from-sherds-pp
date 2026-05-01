import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Plot GA Convergence from CSV.')
parser.add_argument('--csv', type=str, default='result_paper/ga_convergence.csv', help='Path to the input CSV file')
parser.add_argument('--outdir', type=str, default='result_paper', help='Directory to save output plots')
args = parser.parse_args()

outdir = args.outdir
if not os.path.exists(outdir):
    os.makedirs(outdir)

csv_path = args.csv
if not os.path.exists(csv_path):
    print(f"Error: {csv_path} not found.")
    exit(1)

df = pd.read_csv(csv_path)

# ── 1. Overall Convergence ──────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
fig.suptitle('Genetic Algorithm Convergence', fontsize=16, fontweight='bold')

axes[0].plot(df['generation'], df['best_fitness'], label='Best Fitness', linewidth=2, color='#1565C0')
axes[0].set_ylabel('Score', fontsize=12)
axes[0].grid(True, linestyle='--', alpha=0.4)
axes[0].legend(loc='lower right', fontsize=11)

axes[1].plot(df['generation'], df['avg_fitness'], label='Average Fitness', linewidth=2, color='#43A047', linestyle='--')
axes[1].set_ylabel('Score', fontsize=12)
axes[1].grid(True, linestyle='--', alpha=0.4)
axes[1].legend(loc='lower right', fontsize=11)

axes[2].plot(df['generation'], df['worst_fitness'], label='Worst Fitness', linewidth=2, color='#E53935', linestyle=':')
axes[2].set_xlabel('Generation', fontsize=12)
axes[2].set_ylabel('Score', fontsize=12)
axes[2].grid(True, linestyle='--', alpha=0.4)
axes[2].legend(loc='lower right', fontsize=11)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'fitness_curves.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fitness_curves.png")

# ── 2. Individual Component Subplots ────────────────────────────────
components = [
    ('inlier_reward',            'Inlier Reward',             '#1565C0', False),
    ('connectivity_reward',      'Connectivity Reward',       '#2E7D32', False),
    ('cycle_penalty',            'Cycle Penalty',             '#E65100', True),
    ('edge_residual_penalty',    'Edge Residual Penalty',     '#AD1457', True),
    ('rot_residual_penalty',     'Rotation Residual Penalty', '#6A1B9A', True),
    ('overlap_penalty',          'Overlap Penalty',           '#C62828', True),
]

fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
fig.suptitle('Fitness Component Breakdown (Best Individual)', fontsize=16, fontweight='bold', y=0.98)

for ax, (col, title, color, is_penalty) in zip(axes.flat, components):
    values = df[col]
    ax.plot(df['generation'], values, linewidth=2, color=color)
    ax.fill_between(df['generation'], 0, values, alpha=0.15, color=color)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_ylabel('Value', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.axhline(y=0, color='gray', linewidth=0.5, linestyle='-')

    # Annotate final value
    final_val = values.iloc[-1]
    ax.annotate(f'{final_val:.1f}', xy=(df['generation'].iloc[-1], final_val),
                fontsize=9, color=color, fontweight='bold',
                ha='right', va='bottom' if final_val >= 0 else 'top')

axes[-1, 0].set_xlabel('Generation', fontsize=11)
axes[-1, 1].set_xlabel('Generation', fontsize=11)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(outdir, 'fitness_breakdown.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fitness_breakdown.png")

# ── 3. Structural Metrics ──────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
fig.suptitle('Structural Metrics (Best Individual)', fontsize=14, fontweight='bold')

axes[0].plot(df['generation'], df['active_pairs'], linewidth=2, color='#0277BD', marker='o', markersize=3)
axes[0].set_title('Active Pairs', fontsize=12)
axes[0].set_xlabel('Generation')
axes[0].set_ylabel('Count')
axes[0].grid(True, linestyle='--', alpha=0.4)

axes[1].plot(df['generation'], df['largest_component'], linewidth=2, color='#00695C', marker='s', markersize=3)
axes[1].set_title('Largest Connected Component', fontsize=12)
axes[1].set_xlabel('Generation')
axes[1].set_ylabel('Size')
axes[1].grid(True, linestyle='--', alpha=0.4)

axes[2].plot(df['generation'], df['num_components'], linewidth=2, color='#BF360C', marker='^', markersize=3)
axes[2].set_title('Number of Components', fontsize=12)
axes[2].set_xlabel('Generation')
axes[2].set_ylabel('Count')
axes[2].grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(outdir, 'structural_metrics.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved structural_metrics.png")

print(f"\nAll plots saved to {outdir}/")
