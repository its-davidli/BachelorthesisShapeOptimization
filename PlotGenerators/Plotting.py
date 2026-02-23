import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.legend_handler import HandlerTuple
data_circle = np.loadtxt("FinalResults/RectangleToCircle/Figures_and_Data/objective_values.txt", skiprows=1)
data_2Sides = np.loadtxt("FinalResults/Rectangle2Sides/Figures_and_Data/objective_values.txt", skiprows=1)

gradient_norms_circle = np.sqrt(data_circle[:, 2])
gradient_norms_2Sides = np.sqrt(data_2Sides[:, 2])
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
p1, = ax1.plot(scaled_iterations_circle, objective_values_circle, label="Objective Values (Hedgehog Defect)", color=color1)
p2, = ax1.plot(scaled_iterations_2Sides, objective_values_2Sides, label="Objective Values (Uniform Directorfield)", color=color1, linestyle='dashed')
# Plot with log scales on both axes
ax1.set_yscale('log')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xticks(scaled_iterations_circle)
ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel("Shape Gradient Norm Squared", color=color2)
ax2.set_yscale('log')
p3, = ax2.plot(scaled_iterations_circle, gradient_norms_circle, label="Gradient Norms (Hedgehog Defect)", color=color2)
p4, = ax2.plot(scaled_iterations_2Sides, gradient_norms_2Sides, label="Gradient Norms (Uniform Directorfield)", color=color2, linestyle='dashed')
ax2.tick_params(axis='y', labelcolor=color2)
ax2.set_xticks(scaled_iterations_2Sides)
ax2.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
# fig.legend([(p1, p3), (p2,p4)], ['Hedgehog Defect', 'Uniform Directorfield'],
            #    handler_map={tuple: HandlerTuple(ndivide=None)})
fig.legend(fontsize='x-small', loc='upper right', bbox_to_anchor=(0.85, 0.95))
plt.xticks(np.arange(0, maxIterations +1,maxIterations/5), labels= np.arange(0, 101,20))
fig.tight_layout()
    
plt.savefig("PlotGenerators/ObjectiveValue_vs_Iteration.png")
