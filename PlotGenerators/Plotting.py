import matplotlib.pyplot as plt 
import numpy as np

data_circle = np.loadtxt("PreliminaryResults/RectangleToCircleNew/Figures_and_Data/objective_values.txt", skiprows=1)
data_2Sides = np.loadtxt("PreliminaryResults/Rectangle2Sides/Figures_and_Data/objective_values.txt", skiprows=1)

rel_changes_circle = data_circle[:, -2]
rel_changes_2Sides = data_2Sides[:, -2]
objective_values_circle = data_circle[:, 1]
objective_values_2Sides = data_2Sides[:, 1]
iterations_circle = data_circle[:, 0]
iterations_2Sides = data_2Sides[:, 0]

print(iterations_2Sides)

# They have different iterations, to plot them together, we can scale it to 100 iterations and plot them together
# On the same figure plot the relative changes in objective value with scale on the lefthand side

maxIterations = np.lcm(len(iterations_circle), len(iterations_2Sides))
scaled_iterations_circle = np.linspace(0, maxIterations, len(iterations_circle))
scaled_iterations_2Sides = np.linspace(0, maxIterations, len(iterations_2Sides))

plt.figure(figsize=(10, 6))
fig, ax1 = plt.subplots()

color1 = 'tab:blue'
ax1.set_xlabel("Iteration")
ax1.set_ylabel("Objective Value", color=color1)
ax1.plot(scaled_iterations_circle, objective_values_circle, label="Rectangle to Circle", color=color1)
ax1.plot(scaled_iterations_2Sides, objective_values_2Sides, label="Rectangle to 2 Sides", color=color1, linestyle='dashed')
# Plot with log scales on both axes
ax1.set_yscale('log')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xticks(scaled_iterations_circle)
ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel("Shape Gradient Norm Squared", color=color2)
ax2.set_yscale('log')
ax2.plot(scaled_iterations_circle, rel_changes_circle, label="Rectangle to Circle", color=color2)
ax2.plot(scaled_iterations_2Sides, rel_changes_2Sides, label="Rectangle to 2 Sides", color=color2, linestyle='dashed')
ax2.tick_params(axis='y', labelcolor=color2)
ax2.set_xticks(scaled_iterations_2Sides)
ax2.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

fig.tight_layout()
    
plt.savefig("PlotGenerators/ObjectiveValue_vs_Iteration.png")
