import matplotlib.pyplot as plt
import numpy as np
from fenics import *
import dolfin
import numpy as np
from ufl import nabla_div, VectorElement, FiniteElement, MixedElement, split
import math 
import sys
import os
import xml.etree.ElementTree as ET

#Create the new mesh
mesh = Mesh()

# #Read Mesh into the empty mesh
with XDMFFile(MPI.comm_world, "Geometries/OneSidedWallBox/OneSidedWallBox_3D.xdmf") as xdmf_infile:
    xdmf_infile.read(mesh)


def read_xdmf_time_values(xdmf_filename):
    """
    Reads time values from an XDMF file and returns a dictionary
    mapping grid names (e.g., "Vm_0", "Vm_1") to their corresponding time values.
    From https://fenicsproject.discourse.group/t/read-checkpoint-retrieve-checkpoint-time/14384/4
    """
    time_values = []
    
    # Parse the XML tree
    tree = ET.parse(xdmf_filename)
    root = tree.getroot()
    
    # Find all grids inside the temporal collection
    
    for grid in root.findall(".//Grid[@GridType='Uniform']"):
        time_element = grid.find("Time")
        
        if time_element is not None:
            time_values.append(float(time_element.get("Value")))
    
    return time_values

time_values = read_xdmf_time_values("PreliminaryResults/Test121/directors.xdmf")
time_values = time_values[:-1]
for time in time_values:
    print(time)
    S = VectorFunctionSpace(mesh, "CG", 1)
    u = Function(S)

    with XDMFFile(dolfin.MPI.comm_world, "PreliminaryResults/Test121/displacements.xdmf") as xdmf_outfile:
        xdmf_outfile.read_checkpoint(u,"displacement",int(time))

    ALE.move(mesh, u)

director = Function(S)

with XDMFFile(dolfin.MPI.comm_world, "PreliminaryResults/Test121/directors.xdmf") as xdmf_outfile:
    xdmf_outfile.read_checkpoint(director,"director",int(time_values[-1]))

# Now we read the final director in. Now I want to plot the director vectors along the midline of the box. I can do this by creating a line along the midline and then evaluating the director at points along this line.

# Iterate through the vertices of the mesh and find those that are close to the midline (x=0.5, y=0.5). 
midline_points = np.array([]).reshape(0, 3)
midline_director = np.array([]).reshape(0, 3)
vertex_values_u = director.compute_vertex_values()
num_dofs_per_component = int(S.dim()/S.num_sub_spaces())
num_sub_spaces = S.num_sub_spaces()
vector = np.zeros((num_sub_spaces, num_dofs_per_component))
for i in range(num_sub_spaces):
    vector[i] = director.sub(i, deepcopy=True).vector().get_local()
x = S.sub(0).collapse().tabulate_dof_coordinates()
vector = vector.T
r = 0.05
for coord, vec in zip(x, vector):
    if abs(coord[0] - 0.5) < r and abs(coord[1] - 0.5) < r:
        midline_points = np.append(midline_points, np.array([coord]), axis=0)
        midline_director = np.append(midline_director, np.array([vec/np.linalg.norm(vec)]), axis=0)
        

# Sort the midline points and directors by the z-coordinate
sorted_indices = np.argsort(midline_points[:, 2])
midline_points = midline_points[sorted_indices]
midline_director = midline_director[sorted_indices]


# The director has a pi symmetry, so we need to ensure that all vectors point in the same general direction for better visualization. We can do this by checking the dot product with the previous vector.
for i in range(1, midline_director.shape[0]):
    if np.dot(midline_director[i], midline_director[i-1]) < 0:
        midline_director[i] = -midline_director[i]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(midline_points[:,0], midline_points[:,1], midline_points[:,2], midline_director[:,0]/5, midline_director[:,1]/5,   midline_director[:,2]/5)
# ax.plot(midline_points[:,0] + midline_director[:,0], midline_points[:,1] + midline_director[:,1], midline_points[:,2] + midline_director[:,2], 'r-')
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([0, 5])
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.title('Director Field along Midline of the Box')
plt.savefig('director_midline.png')

# calculate the azimuthal and polar angles of the director vectors along the midline
azimuthal_angles = np.arctan2(midline_director[:, 1], midline_director[:, 0])

polar_angles = np.arccos(midline_director[:, 2])
# modulo the pi symmetry, we can take the polar angle to be in the range [0, pi/2]
azimuthal_angles = np.mod(azimuthal_angles, np.pi)

plt.figure()
plt.plot(midline_points[:, 2], azimuthal_angles, label='Azimuthal Angle', color = 'blue')
plt.plot(midline_points[:, 2], (np.pi / 4) * 0.2 * midline_points[:, 2], label='Target Azimuthal Angle', linestyle='--', color = 'blue')
plt.plot(midline_points[:, 2], polar_angles, label='Polar Angle', color = 'red')
# mark pi/4 at the y axis with an extra tick and label it pi/4
plt.yticks(np.arange(0, np.pi/2 + np.pi/8, np.pi/8).tolist() , labels=[r'$0^\circ$', r'$22.5^\circ$', r'$45^\circ$', r'$67.5^\circ$', r'$90^\circ$'])

plt.xlabel('Z coordinate')
plt.ylabel('Angle')
plt.title('Azimuthal and Polar Angles of Director along z-Midline')
plt.legend()
plt.savefig('director_angles_midline.png')