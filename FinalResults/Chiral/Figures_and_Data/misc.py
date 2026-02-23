from fenics import *
from dolfin import *
from ufl_legacy import nabla_div, nabla_grad, VectorElement, FiniteElement, MixedElement, split, atan_2, replace, Jacobian,  tr, variable, shape, Max, sqrt, Min, det, Identity, max_value, min_value
import numpy as np
import math
import matplotlib.pyplot as plt

def print_objective_and_gradient_terms(terms):
    print("\n{:35s} {:>20s} {:>20s} {:>15s} {:>25s} {:>25s}".format("Term Name", "Objective Value", "dJdS Value", "Coeff", "Coeff*Objective", "Coeff*dJdS"))
    print("-"*145)
    for name, obj_expr, grad_expr, coeff in terms:
        obj_val = assemble(obj_expr)
        grad_val = assemble(grad_expr).norm('l2')
        print("{:35s} {:20.6e} {:20.6e} {:15.6e} {:25.6e} {:25.6e}".format(name, obj_val, grad_val, coeff, obj_val*coeff, coeff*grad_val))

def plot_boundary_subdomain(mesh, boundaries, surf_markers_anchoring, save_dir):
    # Visualize boundary subdomain with subdomain_ids=surf_markers, for visual verification of anchoring or moving boundary
    # 0 corresponds to non-active boundary facets, 1 to active boundary facets
    boundary_marked = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundary_marked.set_all(0)
    for marker in surf_markers_anchoring:
        for facet in facets(mesh):
            if boundaries[facet.index()] == marker:
                boundary_marked[facet.index()] = 1
    File(save_dir) << boundary_marked

def write_objective_terms_to_file(save_dir, dic):
    # dictionary dic with keyname in the first row and the value the array containing the values
    # Write dictionary values to file with keys as header
    keys = list(dic.keys())
    num_rows = len(next(iter(dic.values())))
    with open(save_dir + "/Figures_and_Data/objective_values.txt", "w") as f:
        f.write("Iteration\t" + "\t".join(keys) + "\n")
        for i in range(num_rows):
            row = [str(dic[key][i]) for key in keys]
            f.write(f"{i}\t" + "\t".join(row) + "\n")

def plotResults(save_dir, objective_values, shape_gradient_norms, rel_changes, abs_changes, **kwargs):
    
    iterations = np.arange(len(objective_values))

    fig, ax1 = plt.subplots()

    color1 = 'tab:blue'
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Objective Value", color=color1)
    ax1.plot(iterations, objective_values, label="Objective Value", color=color1)
    # Plot with log scales on both axes
    ax1.set_yscale('log')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.set_xticks(iterations)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.set_ylabel("Shape Gradient Norm Squared", color=color2)
    ax2.plot(iterations, shape_gradient_norms, label="Shape Gradient Norm Squared", color=color2)
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_xticks(iterations)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    fig.tight_layout()
    plt.savefig(save_dir + "/Figures_and_Data/objective_and_gradient_norms_log.png")
    plt.close()

    # Plot with linear scales on both axes
    fig_lin, ax1_lin = plt.subplots()
    color1_lin = 'tab:blue'
    ax1_lin.set_xlabel("Iteration")
    ax1_lin.set_ylabel("Objective Value", color=color1_lin)
    ax1_lin.plot(iterations, objective_values, label="Total Objective Value", color=color1_lin)
    ax1_lin.tick_params(axis='y', labelcolor=color1_lin)
    ax1_lin.set_xticks(iterations)
    ax1_lin.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax2_lin = ax1_lin.twinx()
    color2_lin = 'tab:red'
    ax2_lin.set_ylabel("Shape Gradient Norm Squared", color=color2_lin)
    ax2_lin.plot(iterations, shape_gradient_norms, label="Shape Gradient Norm Squared", color=color2_lin)
    ax2_lin.tick_params(axis='y', labelcolor=color2_lin)
    ax2_lin.set_xticks(iterations)
    ax2_lin.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    fig_lin.tight_layout()
    plt.savefig(save_dir + "/Figures_and_Data/objective_and_gradient_norms_linear.png")
    plt.close()

    # plot rel_changes and abs_changes
    fig_lin, ax1_lin = plt.subplots()
    color1_lin = 'tab:blue'
    ax1_lin.set_xlabel("Iteration")
    ax1_lin.set_ylabel("Relative Change", color=color1_lin)
    ax1_lin.plot(iterations, rel_changes, label="Relative Change in Objective", color=color1_lin)
    ax1_lin.tick_params(axis='y', labelcolor=color1_lin)
    ax1_lin.set_xticks(iterations)
    ax1_lin.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax2_lin = ax1_lin.twinx()
    color2_lin = 'tab:red'
    ax2_lin.set_ylabel("Absolute Change", color=color2_lin)
    ax2_lin.plot(iterations, abs_changes, label="Absolute Change in Objective", color=color2_lin)
    ax2_lin.tick_params(axis='y', labelcolor=color2_lin)
    ax2_lin.set_xticks(iterations)
    ax2_lin.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    fig_lin.tight_layout()
    plt.savefig(save_dir + "/Figures_and_Data/rel_changes.png")
    plt.close()

    # plot rel_changes and abs_changes with log scale
    fig_log, ax1_log = plt.subplots()
    color1_log = 'tab:blue'
    ax1_log.set_xlabel("Iteration")
    ax1_log.set_ylabel("Relative Change", color=color1_log)
    ax1_log.plot(iterations, rel_changes, label="Relative Change in Objective", color=color1_log)
    ax1_log.set_yscale('log')
    ax1_log.tick_params(axis='y', labelcolor=color1_log)
    ax1_log.set_xticks(iterations)
    ax1_log.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax2_log = ax1_log.twinx()
    color2_log = 'tab:red'
    ax2_log.set_ylabel("Absolute Change", color=color2_log)
    ax2_log.plot(iterations, abs_changes, label="Absolute Change in Objective", color=color2_log)
    ax2_log.set_yscale('log')
    ax2_log.tick_params(axis='y', labelcolor=color2_log)
    ax2_log.set_xticks(iterations)
    ax2_log.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    
    fig_log.tight_layout()
    plt.savefig(save_dir + "/Figures_and_Data/rel_changes_log.png")
    plt.close()

    for key, values in kwargs.items():
        plt.figure()
        plt.plot(iterations, values, label=key)
        plt.xlabel("Iteration")
        plt.legend()
        plt.xticks(iterations)
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.savefig(save_dir + f"/Figures_and_Data/{key}.png")
        plt.close()

        # plot with log scale
        plt.figure()
        plt.plot(iterations, values, label=key)
        plt.yscale('log')
        plt.xlabel("Iteration")
        plt.legend()
        plt.xticks(iterations)
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.savefig(save_dir + f"/Figures_and_Data/{key}_log.png")
        plt.close()


def plotGeometricalInformation(save_dir, radii, variances_radius, center_of_masses):
    iterations = np.arange(len(radii))

    # Plot radii
    plt.figure()
    plt.plot(iterations, radii, label="Radii")
    plt.xlabel("Iteration")
    plt.legend()
    plt.xticks(iterations)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.savefig(save_dir + "/Figures_and_Data/radii.png")
    plt.close()

    # Plot variances of radii
    plt.figure()
    # plt.plot(iterations, variances_radius, label="Variances of Radii")
    plt.plot(iterations, np.sqrt(variances_radius/(np.square(radii))), label="Relative Standard Deviation of Radii")
    plt.xlabel("Iteration")
    plt.ylabel("$\\sigma/\\mathrm{radius}$")
    plt.legend()
    plt.xticks(iterations)
    plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.savefig(save_dir + "/Figures_and_Data/var_of_radii.png")
    plt.close()

    
    # Plot Center of masses
    center_of_masses = np.array(center_of_masses)
    if np.shape(center_of_masses)[1] == 3:
        x, y, z = center_of_masses[:,0], center_of_masses[:,1], center_of_masses[:,2]
        ax = plt.figure().add_subplot(projection='3d')
        ax.plot(x, y, z, marker='o')
        # Plot arrows
        u = np.diff(x)
        v = np.diff(y)
        w = np.diff(z)
        pos_x = x[:-1] + u/2
        pos_y = y[:-1] + v/2
        pos_z = z[:-1] + w/2

        norm = np.sqrt(u**2+v**2 + w**2) 

        # ax.quiver(pos_x, pos_y, pos_z, u/norm, v/norm, w/norm, zorder=5, pivot="middle")
        plt.savefig(save_dir + "/Figures_and_Data/center_of_masses.png")
        plt.close()

    if np.shape(center_of_masses)[1] == 2:
        x, y= center_of_masses[:,0], center_of_masses[:,1]
        fig, ax = plt.subplots()
        ax.plot(x, y, marker='o')
        # Plot arrows
        u = np.diff(x)
        v = np.diff(y)
        pos_x = x[:-1] + u/2
        pos_y = y[:-1] + v/2

        norm = np.sqrt(u**2+v**2) 

        ax.quiver(pos_x, pos_y, u/norm, v/norm, zorder=5, pivot="middle")
        plt.savefig(save_dir + "/Figures_and_Data/center_of_masses.png")
        plt.close()

def center_of_mass(mesh, d, dx):
    x = Expression('x[0]', degree = 1)
    y = Expression('x[1]', degree = 1)
    x_coord = assemble(x * dx) / assemble(1.0 * dx)
    y_coord = assemble(y * dx) / assemble(1.0 * dx)
    center_of_mass = [x_coord, y_coord]
    if d == 3: 
        z = Expression('x[2]', degree = 1)
        z_coord = assemble(z * dx) / assemble(1.0 * dx)
        center_of_mass = [x_coord, y_coord, z_coord]
    return center_of_mass


def norm_variance_radius(mesh, surf_markers, boundaries,d, dx):
    com = center_of_mass(mesh, d, dx)
    # Compute variance of the radius from the center of mass
    # Only consider boundary points for variance calculation
    boundary_indices = []
    for facet in facets(mesh):
        if boundaries[facet.index()] in surf_markers:
            for vertex in vertices(facet):
                boundary_indices.append(vertex.index())
    boundary_indices = np.unique(boundary_indices)
    coords = mesh.coordinates()
    boundary_coords = coords[boundary_indices]
    radii = np.linalg.norm(boundary_coords - com, axis=1)
    return np.mean(radii), np.var(radii)