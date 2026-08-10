import pandas as pd
import matplotlib.pyplot as plt

# Load data
df_a = pd.read_csv('POT_A/result_paper/ga_convergence.csv')
df_e = pd.read_csv('POT_E/result_paper/ga_convergence.csv')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Pot A
axes[0].plot(df_a['generation'], df_a['best_fitness'], label='Best Fitness', color='blue', linewidth=2)
axes[0].plot(df_a['generation'], df_a['avg_fitness'], label='Avg Fitness', color='orange', alpha=0.7)
axes[0].set_title('Pot A Convergence (< 10 gens)')
axes[0].set_xlabel('Generation')
axes[0].set_ylabel('Fitness')
axes[0].grid(True, linestyle='--', alpha=0.6)
axes[0].legend()

# Pot E
axes[1].plot(df_e['generation'], df_e['best_fitness'], label='Best Fitness', color='red', linewidth=2)
axes[1].plot(df_e['generation'], df_e['avg_fitness'], label='Avg Fitness', color='orange', alpha=0.7)
axes[1].set_title('Pot E Convergence (Plateaus)')
axes[1].set_xlabel('Generation')
axes[1].set_ylabel('Fitness')
axes[1].grid(True, linestyle='--', alpha=0.6)
axes[1].legend()

plt.tight_layout()
plt.savefig('ieee_paper_submission/convergence.png', dpi=300)
print("Plot saved to ieee_paper_submission/convergence.png")
