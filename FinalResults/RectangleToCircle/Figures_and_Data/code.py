from fenics import *
from dolfin import *
import numpy as np
from ufl_legacy import nabla_div, nabla_grad, VectorElement, FiniteElement, MixedElement, split, atan_2, replace, Jacobian
import math 
import sys
import os
from IPython import embed
import matplotlib.pyplot as plt
import scipy
import yaml # For reading configuration files
from DeliverableShapeOptimizationLCsMethods import *
from DeliverableShapeOptimizationMisc import *
# set_log_level(LogLevel.WARNING)
    
# Load configuration from YAML file
with open("config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

# Code to compute the director field n as a minimizer of the Landau de Gennes (LdG) energy functional
# using the discontinuous Galerkin finite element method. The code takes as input a Fenics/dolfin suitable geometry
# in xdmf format and ratio of constants L/c and a/c, while assuming b = c.

#Create folder to save results
save_fold = config['save_fold']
save_subdir = sys.argv[1] 
save_dir = "./"+save_fold + "/" + save_subdir
os.system("mkdir "+" "+save_dir)
os.system("mkdir "+" "+save_dir + "/Figures_and_Data")

# Save the configuration file and code to the save directory
os.system("cp config.yaml " + save_dir + "/Figures_and_Data/config.yaml")
os.system("cp "+ __file__ + " " + save_dir + "/Figures_and_Data/code.py")
os.system("cp "+ 'DeliverableShapeOptimizationLCsMethods.py' + " " + save_dir + "/Figures_and_Data/methods.py")
os.system("cp "+ 'DeliverableShapeOptimizationMisc.py' + " " + save_dir + "/Figures_and_Data/misc.py")

# Folder for initial mesh files
subfolder = config['subfolder']
geom_folder = config['geom_folder']
mesh_name = config['mesh_name']

# Clear (earlie) debugging file
open(save_dir + '/Figures_and_Data/debugging.txt', 'w').close

# Read finite element parameters from config
finite_element = config['finite_element']
finite_element_degree = config['finite_element_degree']


#Import .xdmf geometry file into mesh + Mesh value collection for boundary conditions
surf_markers_anchoring = tuple(config['surf_markers_anchoring'])
surf_markers_moving = tuple(config['surf_markers_moving'])
mesh = Mesh()

mvc_subdomain = MeshValueCollection("size_t", mesh, mesh.topology().dim())
mvc_boundaries = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)
if config['dimensions'] == 2:
    with XDMFFile(MPI.comm_world, geom_folder+subfolder+"/"+mesh_name+"_2D.xdmf") as xdmf_infile:
        xdmf_infile.read(mesh)
        xdmf_infile.read(mvc_subdomain, "")

    with XDMFFile(MPI.comm_world, geom_folder+subfolder+"/"+mesh_name+"_1D.xdmf") as xdmf_infile:
        xdmf_infile.read(mvc_boundaries, "")

elif config['dimensions'] == 3:
    with XDMFFile(MPI.comm_world, geom_folder+subfolder+"/"+mesh_name+"_3D.xdmf") as xdmf_infile:
        xdmf_infile.read(mesh)
        xdmf_infile.read(mvc_subdomain, "")

    with XDMFFile(MPI.comm_world, geom_folder+subfolder+"/"+mesh_name+"_2D.xdmf") as xdmf_infile:
        xdmf_infile.read(mvc_boundaries, "")

else:
    raise ValueError("Configuration dimensions must be 2 or 3.")

domains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc_boundaries)

if mesh.topology().dim() != config['dimensions']:
    raise ValueError("Mesh dimension does not match configuration dimensions.")
d = mesh.topology().dim()
q_degree = 5 # Quadrature degree for integration
dx = Measure('dx', domain=mesh, subdomain_data=domains) # Standard volume measure with subdomains
dx = dx(metadata={'quadrature_degree': q_degree}) 
ds_controlvariable = Measure('ds', domain=mesh, subdomain_data=boundaries, subdomain_id=surf_markers_moving) # Measure for boundary integrals over boundary that can be deformed by the shape optimization
ds_anchoring = Measure('ds', domain=mesh, subdomain_data=boundaries, subdomain_id=surf_markers_anchoring) # Measure for boundary integrals over boundary with anchoring conditions for LC director field

plot_boundary_subdomain(mesh, boundaries, surf_markers_anchoring, save_dir + "/boundary_subdomain_anchoring.pvd")
plot_boundary_subdomain(mesh, boundaries, surf_markers_moving, save_dir + "/boundary_subdomain_moving.pvd")

# Read physical parameters from config
n_indep_comp = int(d*(d+1)/2 - 1) #Number of independent components for a symmetric traceless tensor)
L_c = float(config["L_c"])
a_B = float(config["a_B"])
if d==2: S0 = np.sqrt(2.0*a_B) # Corresponding LdG S eigenvalue
if d==3: S0 = 1/4*(np.sqrt(24.0*a_B + 1)+1) # Corresponding LdG S eigenvalue

# Set default algorithmic parameters (Armijo line search)
maxIter = int(config['maxIter'])        # Maximum number of shape gradient steps.
sigma = float(config['sigma'])          # Armijo line search slope parameter
beta = float(config['beta'])            # Armijo backtracking parameter
alphaInit = float(config['alphaInit'])  # Initial step size for Armijo line search
alphaMin = float(config['alphaMin'])    # Fallback minimal step size to terminate Armijo line search.
tol = 1E-12                             # General tolerance for floating point comparisons

k_bc = Constant(float(config['k_bc']))  # Penalty constant for the (weak) anchoring term in the LdG energy functional 

# Parameters for the forward problem
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
            "eliminate_zeros": True, \
            "precompute_basis_const": True, \
            "precompute_ip_const": True}

# Vector function space for the mesh deformation
S = VectorFunctionSpace(mesh, "CG", 1)
s = Function(S)

# Define the function space for the state variable q
Lagrange = FiniteElement(finite_element, mesh.ufl_cell(), finite_element_degree) 
mixed_element = MixedElement([Lagrange for i in range(n_indep_comp)]) #Mixed element for q coefficients     
W = FunctionSpace(mesh,mixed_element) 

q_ = Function(W) 
p_ = Function(W)

## Define forward problem 
# Specify Q tensor at boundary
normals = FacetNormal(mesh)

if config['anchoring'] == 'vertical':
    Q_b = S0*(outer(normals,normals)-(1/d)*Identity(d)) 
elif config['anchoring'] == 'planar' and d == 2:
    tangentials = as_vector([-normals[1], normals[0]])
    Q_b = (S0/2.0)*(d*outer(tangentials, tangentials)-Identity(d)) 
else:
    raise NotImplementedError("Anchoring type not implemented or not compatible with dimension.")


Energy = (LG_energy(q_,a_B,L_c, d))*dx + LG_boundaryPenalty(q_,Q_b,k_bc, d) * ds_anchoring
Dv_Energy= derivative(Energy, q_, p_)

# Define the target state variable q_target based on the configuration. This is the state variable that the shape optimization will try to achieve by deforming the domain.
q_target = get_objective(config['target_geometry'], mesh, S0, d)
q_target_proj = project(q_target, W) # Project target into finite element space for visualization


if d == 2:
    File(save_dir + '/target.pvd') << q_target_proj # Save the target Q
else: 
    xdmffile_target = XDMFFile(save_dir + '/target.xdmf')
    xdmffile_target.write(q_target_proj,0)

target_orientation, target_S = compute_orientation(q_target, mesh, d)
File(save_dir + '/target_director.pvd') << target_orientation # Save the target director field
File(save_dir + '/target_S.pvd') << target_S # Save the target scalar order parameter field

# Define the objective function

# The objective function is the sum of the squared difference between the state variable q_ and the target q_target
if d == 2: objective_main = 2*(dot(q_ - q_target, q_ - q_target))/(assemble(1*dx))*dx
if d == 3: objective_main = ((dot(q_ - q_target, q_ - q_target)) + (q_[1]-q_target[1])*(q_[1]-q_target[1]) + (q_[2]-q_target[2])*(q_[2]-q_target[2]) + (q_[4]-q_target[4])*(q_[4]-q_target[4]) + (q_[0]+ q_[3]-q_target[0]-q_target[3])*(q_[0]+ q_[3]-q_target[0]-q_target[3]))/(assemble(1*dx))*dx



# Combine objective with the LdG energy to form the Lagrangian of the shape optimization problem and set up adjoint and forward problems
objective = objective_main
Lagrangian = objective + Dv_Energy
forwardPDE = derivative(Lagrangian, p_, TestFunction(W)) 
adjointPDE = replace(derivative(Lagrangian, q_, TestFunction(W)), {p_: TrialFunction(W)})
forwardJacobian = derivative(forwardPDE, q_)
adjointJacobian = derivative(adjointPDE, p_)


# Compute and save initial guess for the state variable q
initial_guess = compute_initial_guess(mesh, W,config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring'])
assign(q_, initial_guess)
assign(p_, initial_guess)

# Define functions which are going to hold the solution of the shape gradient,
# as well as a scaled version of the negative shape gradient.
shapeGradient = Function(S)
displacement = Function(S)

# Setup the inner product for the shape gradient computation
inner_product = config['inner_product']
if inner_product == "H1":
    displacementInnerProduct = inner(grad(TrialFunction(S)), grad(TestFunction(S))) * dx  
elif inner_product == "elasticity" or inner_product == "elasticity_trace_free":
    strain, C = get_elasticity_operators(config['E'], config['nu'], d, inner_product)
    displacementInnerProduct = inner(C(strain(TrialFunction(S))),strain(TestFunction(S))) * dx 
else: 
    raise ValueError("Inner product not supported")

delta_id = Constant(config['delta']) # Regularization parameter for the identity part of the inner product
delta_beltrami = Constant(config['delta_beltrami']) # Regularization parameter for the Beltrami part of the inner product (if tangential smoothing is used)
displacementInnerProduct += delta_id*inner(TrialFunction(S), TestFunction(S)) * dx
if config['tangential_smoothing']:
    displacementInnerProduct += delta_beltrami*inner(tang_grad(TrialFunction(S), normals), tang_grad(TestFunction(S), normals)) * ds_controlvariable

objective_values, alphas, shape_gradient_norms, rel_changes, abs_changes, objectives_main, objectives_meshquality, volumes, variances_radius, radii, center_of_masses = [], [], [], [], [], [], [], [], [], [] , []


# Set up XDMF files for saving results in 3D
if d == 3:
    xdmffile_results = XDMFFile(save_dir + '/results.xdmf')
    xdmffile_adjoints = XDMFFile(save_dir + '/adjoints.xdmf')
    xdmffile_displacements = XDMFFile(save_dir + '/displacements.xdmf')
    xdmffile_directors = XDMFFile(save_dir + '/directors.xdmf')
# Run a shape gradient loop.
iteration = 0
while iteration < maxIter:
    X = SpatialCoordinate(mesh)

    # Solve the forward PDE.
    # compute initial guess
    if iteration != 0:
        initial_guess =  compute_initial_guess(mesh, W, config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring'])
        assign(q_, initial_guess)
        assign(p_, initial_guess)

    # Solve the forward PDE with relaxation paramaters in the Newton solver
    solveMultRelaxation(config['forward_solver_relaxations'], forwardPDE,0, q_, None, forwardJacobian, ffc_options)
    
    J = assemble(objective)
    # Solve the adjoint PDE.
    solve(lhs(adjointPDE) == rhs(adjointPDE), p_, None)

    # Evaluate the shape derivative. 
    nh = compute_normals(mesh)
    dJds_main = derivative(objective + Dv_Energy, X, TestFunction(S))  
    dJds = dJds_main

    # Riesz representation on the interior: find dJds_1 in S s.t.
    # ∫_∂Ω <dJds_1, v> dx = dJds(v) for all v in S
    u_b = TrialFunction(S)
    v_b = TestFunction(S)
    a_b = inner(u_b, v_b) * dx
    L_b = dJds
    A_b = assemble(a_b, keep_diagonal=True)
    b_b = assemble(L_b)
    # Ensure interior dofs get identity to avoid singular system
    A_b.ident_zeros()
    dJds_rep = Function(S)
    solve(A_b, dJds_rep.vector(), b_b)


    # Riesz representation on the boundary: find dJds_1 in S s.t.
    # ∫_∂Ω <dJds_1, v> ds = dJds(v) for all v in S
    u_b = TrialFunction(S)
    v_b = TestFunction(S)
    a_b = inner(u_b, v_b) * ds_controlvariable
    L_b = dJds
    A_b = assemble(a_b, keep_diagonal=True)
    b_b = assemble(L_b)
    # Ensure interior dofs get identity to avoid singular system
    A_b.ident_zeros()
    dJds_1 = Function(S)
    solve(A_b, dJds_1.vector(), b_b)

    # Use only the normal component of this boundary vector for the shape gradient
    N_vec = CGNormal(mesh)
    dJds_1_n = project(inner(dJds_1, N_vec) * N_vec, S)

    if config['use_normal_dJds_projection']:
        solve(displacementInnerProduct == inner(dJds_1_n,TestFunction(S))*dx, shapeGradient)
    else:
        solve(displacementInnerProduct == dJds, shapeGradient)

    
    # Evaluate the squared norm of the shape gradient induced by the (regularized)
    # elasticity inner product.
    normShapeGradient2 = sum(shapeGradient.vector() * shapeGradient.vector())

    # Export results
    if d == 2:
        File(save_dir + '/forwardSol-{0:03d}.pvd'.format(iteration)) << q_
        File(save_dir + '/adjointSol-{0:03d}.pvd'.format(iteration)) << p_    
    else: 
        xdmffile_results.write(q_,iteration)
        xdmffile_adjoints.write(p_,iteration)
    
    File(save_dir + f"/dJds_rep_{iteration}.pvd") << dJds_rep
    File(save_dir + f"/dJds_boundary_{iteration}.pvd") << dJds_1
    File(save_dir + f"/dJds_1_normal_{iteration}.pvd") << dJds_1_n
    File(save_dir + f"/N_vec_{iteration}.pvd") << N_vec

    M, S_param = compute_orientation(q_, mesh, d)
    File(save_dir + '/directorfield-{0:03d}.pvd'.format(iteration)) << M
    File(save_dir + '/scalarorderparameter-{0:03d}.pvd'.format(iteration)) << S_param
    neg_shape_gradient = Function(S)
    neg_shape_gradient.assign(-shapeGradient)
    File(save_dir + f"/neg_shape_gradient_{iteration}.pvd") << neg_shape_gradient
    plot_boundary_subdomain(mesh, boundaries, surf_markers_anchoring, save_dir + f"/boundary_subdomain_anchoring_iter{iteration}.pvd")

    # Printing information
    terms = [
    ("objective_main", objective_main, dJds_main, 1.0),
    # Add more terms as needed
        ]
    print(f'Iteration {iteration}')
    print_objective_and_gradient_terms(terms)
    objective_values.append(J)
    shape_gradient_norms.append(normShapeGradient2)
    objectives_main.append(assemble(objective_main))
    objectives_meshquality.append(assemble((1/(CellVolume(mesh)+Constant(1e-10))**2)/assemble(1*dx)*dx))
    center_of_masses.append(center_of_mass(mesh, d, dx))
    radii.append(norm_variance_radius(mesh, surf_markers_anchoring, boundaries, d, dx)[0])
    variances_radius.append(norm_variance_radius(mesh, surf_markers_anchoring, boundaries, d, dx)[1])
    volumes.append(assemble(1*dx))


    # Store the mesh associated with the current iterate, as well as its objective value.
    referenceMeshCoordinates = mesh.coordinates().copy()
    # Begin Armijo line search.
    lineSearchSuccessful = False
    sub_iteration = 0
    if iteration == 0:
        alpha = alphaInit
    else:
        alpha = min(config['alphaMaxFactor']*alphaInit, alpha / beta)
        
    while (lineSearchSuccessful == False) and (alpha > alphaMin):
        # Assign the mesh displacement vector field.
        mesh.coordinates()[:] = referenceMeshCoordinates
        displacement.assign(- alpha * shapeGradient)

        # Update the mesh by adding the displacement to the reference mesh.
        ALE.move(mesh, displacement)
        File(save_dir + f"/mesh_iter{iteration}_sub{sub_iteration}.pvd") << mesh
        
        # Solve the forward PDE.
        assign(q_, compute_initial_guess(mesh,W, config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring']))
        # Solve the forward PDE with the updated mesh.
        solveMultRelaxation(config['forward_solver_relaxations'], forwardPDE,0, q_, None, forwardJacobian, ffc_options)

        trialJ = assemble(objective)


        # Evaluate the Armijo condition and reduce the step size if necessary.
            
        # Write debugging information to a file
        if (trialJ <= J - sigma * alpha * normShapeGradient2):
            lineSearchSuccessful = True
            alphas.append(alpha)
            # Plotting in the 3D case
            if d == 3:
                xdmffile_displacements.write_checkpoint(displacement,"displacement",iteration,XDMFFile.Encoding.HDF5, True)
                M, S_param = compute_orientation(q_, mesh, d)
                xdmffile_directors.write_checkpoint(M,"director",iteration,XDMFFile.Encoding.HDF5, True)
        
        else:
            alpha *= beta
        # Writing debugging information to a file
        with open(save_dir + '/Figures_and_Data/debugging.txt', 'a') as debug_file:
            debug_file.write(f"It.: {iteration}\t")
            debug_file.write(f"Alpha (step size): {alpha/beta:9.2e}\t")
            debug_file.write(f"Obj. (J): {J:9.2e}\t")
            debug_file.write(f"Trial Obj. (trialJ): {trialJ:9.2e}\t")
            debug_file.write(f"Norm of Shape Gradient squared: {normShapeGradient2:12.2e}\t")
            debug_file.write(f"Armijo Condition: {lineSearchSuccessful}\n")
            debug_file.write("-" * 40 + "\n")

            
        # Increment the sub-iteration counter.
        sub_iteration += 1
    if lineSearchSuccessful == False: 
        print("Line search failed to find a suitable step size.")
        alphas.append(0.0)
        break
    # Occasionally display some information.


    rel_change = abs(trialJ - objective_values[-1]) / (abs(objective_values[-1]) + 1e-12)
    rel_changes.append(rel_change)
    abs_change = abs(trialJ - objective_values[max(iteration-2,0)])/min(3.0, iteration+1)
    abs_changes.append(abs_change)

    if iteration > 0:
        # Save intermediate results
        write_objective_terms_to_file(save_dir, {'objective_values': objective_values, 'shape_gradient_norms_squared': shape_gradient_norms, 'meshquality': objectives_meshquality, 'alphas': alphas, 'rel_changes': rel_changes, 'abs_changes': abs_changes})
        plotResults(save_dir, objective_values, shape_gradient_norms, rel_changes, abs_changes, alphas=alphas, meshqualities=objectives_meshquality, volumes=volumes)
        plotGeometricalInformation(save_dir, radii,variances_radius, center_of_masses)

        # Check stopping criteria (changes in objective functional below threshold or objective functional increased)
        if rel_change < float(config['rel_change_stopping_value']):
            print(f"Stopping: Relative change in objective ({rel_change:.2e}) is below threshold.")
            break
        if abs_change < float(config['abs_change_stopping_value']):
            print(f"Stopping: Absolute change in objective ({abs_change:.2e}) is below threshold.")
            break
        if trialJ> objective_values[-1]:
            print(f"Stopping: Objective functional increased: {trialJ} > {objective_values[-1]}.")
            break

    # Increment the iteration counter.
    iteration += 1