import collections

import matplotlib.pyplot as plt 

data_everything = []
def read_data_to_dict(filepath):
    data_dict = collections.defaultdict(list)
    with open(filepath, 'r') as file:
        headers = file.readline().strip().split()
        for line in file:
            values = [float(x) for x in line.strip().split()]
            for key, value in zip(headers, values):
                data_dict[key].append(value)
    return dict(data_dict)

data_everything = read_data_to_dict('Results/LCEllipse/mesh_quality_terms.txt')
data_no_beltrami = read_data_to_dict('Results/LCEllipseNoBeltrami/mesh_quality_terms.txt')
data_no_beltrami_no_proj = read_data_to_dict('Results/LCEllipseNoBeltramiNoNormalProj/mesh_quality_terms.txt')
data_no_beltrami_no_proj_mesh_qual = read_data_to_dict('Results/LCEllipseNoBeltramiNoNormalProjMeshQual/mesh_quality_terms.txt')

plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], data_everything['objectives_meshquality'], label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], data_no_beltrami['objectives_meshquality'], label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], data_no_beltrami_no_proj_mesh_qual['objectives_meshquality'], label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], data_no_beltrami_no_proj['objectives_meshquality'], label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('Mesh Quality')
plt.legend(fontsize=11.5)
plt.grid(True)
plt.tight_layout()
plt.savefig('mesh_quality_comparison.png')

# Plotting the main objective values
plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], data_everything['objectives_main'], label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], data_no_beltrami['objectives_main'], label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], data_no_beltrami_no_proj_mesh_qual['objectives_main'], label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], data_no_beltrami_no_proj['objectives_main'], label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('(Main) Objective')
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('main_objective_comparison.png')

# Plotting the variance of the radius
plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], data_everything['variance'], label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], data_no_beltrami['variance'], label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], data_no_beltrami_no_proj_mesh_qual['variance'], label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], data_no_beltrami_no_proj['variance'], label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('Variance of Radii')
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('variance_radius_comparison.png')

# Plotting mesh quality normalized by volume
plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], 
         [q/v for q, v in zip(data_everything['objectives_meshquality'], data_everything['volumes'])], 
         label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], 
         [q/v for q, v in zip(data_no_beltrami['objectives_meshquality'], data_no_beltrami['volumes'])], 
         label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], 
         [q/v for q, v in zip(data_no_beltrami_no_proj_mesh_qual['objectives_meshquality'], data_no_beltrami_no_proj_mesh_qual['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], 
         [q/v for q, v in zip(data_no_beltrami_no_proj['objectives_meshquality'], data_no_beltrami_no_proj['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('Mesh Quality / Volume')
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('mesh_quality_comparison_normalized.png')

# Plotting main objective normalized by volume
plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], 
         [o/v for o, v in zip(data_everything['objectives_main'], data_everything['volumes'])], 
         label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], 
         [o/v for o, v in zip(data_no_beltrami['objectives_main'], data_no_beltrami['volumes'])], 
         label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], 
         [o/v for o, v in zip(data_no_beltrami_no_proj_mesh_qual['objectives_main'], data_no_beltrami_no_proj_mesh_qual['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], 
         [o/v for o, v in zip(data_no_beltrami_no_proj['objectives_main'], data_no_beltrami_no_proj['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('(Main) Objective / Volume')
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('main_objective_comparison_normalized.png')

# Plotting variance of radius normalized by volume
plt.figure(figsize=(10, 6))
plt.plot(data_everything['Iteration'], 
         [var/v for var, v in zip(data_everything['variance'], data_everything['volumes'])], 
         label='With Projection onto the Normal and boundary smoothing')
plt.plot(data_no_beltrami['Iteration'], 
         [var/v for var, v in zip(data_no_beltrami['variance'], data_no_beltrami['volumes'])], 
         label='With Projection onto the Normal, no boundary smoothing')
plt.plot(data_no_beltrami_no_proj_mesh_qual['Iteration'], 
         [var/v for var, v in zip(data_no_beltrami_no_proj_mesh_qual['variance'], data_no_beltrami_no_proj_mesh_qual['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing, with mesh quality penalty')
plt.plot(data_no_beltrami_no_proj['Iteration'], 
         [var/v for var, v in zip(data_no_beltrami_no_proj['variance'], data_no_beltrami_no_proj['volumes'])], 
         label='No Projection onto the Normal, no boundary smoothing')
plt.xlabel('Iterations')
plt.ylabel('Variance of Radii / Volume')
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('variance_radius_comparison_normalized.png')