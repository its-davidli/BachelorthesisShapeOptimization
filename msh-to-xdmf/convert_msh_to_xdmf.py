import sys
import meshio
import numpy as np
import os

geom_folder = "Geometries/"
subfolder = sys.argv[2] + "/"
mesh_name = sys.argv[2]
msh = meshio.read("msh-to-xdmf/"+sys.argv[1]+".msh")
os.system("mkdir "+" "+ geom_folder+subfolder)
os.system("mkdir "+" "+ geom_folder+subfolder)
dimension = int(sys.argv[3])
if dimension == 2:
    body_element_tag = "triangle"
    TwoDpoints= msh.points
    TwoDpoints = np.delete(TwoDpoints, obj=2, axis=1) # Removes the 2D points'

    line_data = msh.cell_data_dict["gmsh:physical"][body_element_tag]
    meshio.write(geom_folder+subfolder+"/"+mesh_name+"_2D.xdmf",
        meshio.Mesh(points=TwoDpoints,
        cells={body_element_tag: msh.cells_dict[body_element_tag]},
        cell_data={"bnd_marker": [line_data]})
        ) 

    body_element_tag = "line"



    line_data = msh.cell_data_dict["gmsh:physical"][body_element_tag]
    meshio.write(geom_folder+subfolder+"/"+mesh_name+"_1D.xdmf",
        meshio.Mesh(points=TwoDpoints,
        cells={body_element_tag: msh.cells_dict[body_element_tag]},
        cell_data={"bnd_marker": [line_data]})
        ) 

elif dimension == 3:
    body_element_tag = "tetra"
    TwoDpoints= msh.points
    # TwoDpoints = np.delete(TwoDpoints, obj=2, axis=1) # Removes the 2D points'

    line_data = msh.cell_data_dict["gmsh:physical"][body_element_tag]
    meshio.write(geom_folder+subfolder+"/"+mesh_name+"_3D.xdmf",
        meshio.Mesh(points=TwoDpoints,
        cells={body_element_tag: msh.cells_dict[body_element_tag]},
        cell_data={"bnd_marker": [line_data]})
        ) 

    body_element_tag = "triangle"



    line_data = msh.cell_data_dict["gmsh:physical"][body_element_tag]
    meshio.write(geom_folder+subfolder+"/"+mesh_name+"_2D.xdmf",
        meshio.Mesh(points=TwoDpoints,
        cells={body_element_tag: msh.cells_dict[body_element_tag]},
        cell_data={"bnd_marker": [line_data]})
        ) 