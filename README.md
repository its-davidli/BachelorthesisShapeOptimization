# BachelorthesisShapeOptimization
This is the codebase of the Bachelor Thesis: *Controlling electric and nematic fields through shape optimization.* supervised by  Prof. Dr. Ulrich Schwarz and Santiago Gomez Melo at the Department of Physics and Astronomy, Institute for Theoretical Physics (Heidelberg University).

To reproduce the results, one can find the code file *code.py* in the single results folders, including the config file *config.yaml*. The used Geometries are in the Geometries Folder. The docker image ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30 was used.

For future work, the starting point should be the code DeliverableShapeOptimizationLCs.py with the dependencies DeliverableShapeOptimizationLCsMethods.py
 

*library versions:*  
pip list --format=columns | grep -E "fenics|dolfin|dolfin-adjoint|numpy"

dolfin-adjoint       2023.3.0  
fenics-dijitso       2019.2.0.dev0  
fenics-dolfin        2019.2.0.dev0  
fenics-ffc           2019.2.0.dev0  
fenics-fiat          2019.2.0.dev0  
fenics-ufl-legacy    2022.3.0  
numpy                1.21.5