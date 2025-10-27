from fenics import *
from dolfin import *
from dolfin_adjoint import *
#from create_mesh import c_x, c_y
# import moola
#from scipy.optimize import minimize
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
os.system("mkdir "+" "+save_dir + "/Figures")

# Save the configuration file in the save directory
os.system("cp config.yaml " + save_dir + "/config.yaml")
os.system("cp "+ __file__ + " " + save_dir + "/code.py")
os.system("cp "+ 'DeliverableShapeOptimizationLCsMethods.py' + " " + save_dir + "/methods.py")

# Folder for initial mesh files
subfolder = config['subfolder']
geom_folder = config['geom_folder']
mesh_name = config['mesh_name']

# Clear debugging file
open(save_dir + '/debugging.txt', 'w').close

finite_element = config['finite_element']
finite_element_degree = config['finite_element_degree']
surf_markers = config['surf_markers']

#Import .xdmf geometry file into mesh + Mesh value collection for boundary conditions
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

domains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc_boundaries)

dx = Measure('dx', domain=mesh, subdomain_data=domains)
q_degree = 5
dx = dx(metadata={'quadrature_degree': q_degree}) 
ds_controlvariable = Measure('ds', domain=mesh, subdomain_data=boundaries)
ds_anchoring = Measure('ds', domain=mesh, subdomain_data=boundaries, subdomain_id=surf_markers[0])
d = mesh.topology().dim()
target_geometry = config['target_geometry']

n_indep_comp = int(d*(d+1)/2 - 1) #Number of independent components for a symmetric traceless tensor)
L_c = float(config["L_c"])
a_B = float(config["a_B"])
S0 = np.sqrt(2.0*a_B) # Corresponding LdG S eigenvalue
β = float(config["beta_param"])
x0 = float(config["x0"])

# Set default algorithmic parameters (Armijo line search)
maxIter = int(config['maxIter'])        # Maximum number of shape gradient steps.
sigma = float(config['sigma'])          # Armijo line search slope parameter
beta = float(config['beta'])            # Armijo backtracking parameter
alphaInit = float(config['alphaInit'])  # Initial step size for Armijo line search
alphaMin = float(config['alphaMin'])    # Fallback minimal step size to terminate Armijo line search.

# Penalty constants for constraints and terms in the objective
k_vol = Constant(float(config['k_vol']))                                    # Penalty constant for volume constraint
k_com = Constant(float(config['k_com']))                                    # Penalty constant for center of mass constraint
k_bc = Constant(float(config['k_bc']))                                      # Penalty constant for the anchoring term 
k_boundarylength = Constant(float(config['k_boundarylength']))                          # Penalty constant for the boundary length constraint
k_meshquality = Constant(float(config['k_meshquality']))                    # Penalty constant for the mesh quality penalty
eps_meshquality = Constant(float(config['eps_meshquality']))                # Regularize the meshquality term
k_edge = Constant(float(config['k_edge']))                                  # Penalty constant for the maximum cell edge length constraint

# Define data and auxiliary functions for the elasticity inner product.
E = Constant(float(config['E'])) # Young's modulus
nu = Constant(float(config['nu'])) # Poisson's ratio
lmbda = nu*E/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))
def strain(u): return sym(nabla_grad(u))
def C(epsilon): return lmbda * epsilon + 2 * mu * tr(epsilon) * Identity(d)



x = Expression('x[0]', degree = 1)
y = Expression('x[1]', degree = 1)
if d==3: z = Expression('x[2]', degree = 1)


Vol0 = Constant(assemble(1.0*dx)) # Initial/Target volume of the domain
c_x0 = Constant(assemble(x*dx)/Vol0) # Initial/Target center of mass x coordinate
c_y0 = Constant(assemble(y*dx)/Vol0) # Initial/Target center of mass y coordinate
Boundary0 = Constant(assemble(1.0*ds_controlvariable)) # Initial/Target boundary length of the domain


parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
            "eliminate_zeros": True, \
            "precompute_basis_const": True, \
            "precompute_ip_const": True}


# Vector function space for the mesh deformation
S = VectorFunctionSpace(mesh, "CG", 2)
s = Function(S)

current_volume = Constant(assemble(1.0*dx))
current_boundary_length = Constant(assemble(1.0*ds_controlvariable))

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
    Q_b = (S0/2.0)*(d*outer(normals,normals)-Identity(d)) 
elif config['anchoring'] == 'planar':
    if d == 2:
        tangentials = as_vector([-normals[1], normals[0]])
        Q_b = (S0/2.0)*(d*outer(tangentials, tangentials)-Identity(d)) 
    else:
        raise NotImplementedError("Tangential vector computation is only implemented for 2D meshes.")

Energy = (LG_energy(q_,a_B,L_c,β,x0, d))*dx + LG_boundaryPenalty(q_,Q_b,k_bc, d) * ds_anchoring
Dv_Energy= derivative(Energy, q_, p_)


# Define the target state variable q_target
if d == 2 and target_geometry == "circle":
    theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
    q_target = Expression(('S0*(cos(theta)*cos(theta)-0.5)', 'S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1 = q_target_proj.split()

elif d == 2 and target_geometry == "uniform_vertical":
    q_target = Expression(('-S0*(0.5)', '0'), S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1 = q_target_proj.split()

elif d == 3 and target_geometry == "sphere":
    phi = Expression('atan2((x[1]),(x[0]))', degree = 1)
    theta = Expression('acos(x[2]/std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+1e-20))', degree = 1)
    q_target = Expression(('S0*(cos(phi)*cos(phi)*sin(theta)*sin(theta)-0.5)','S0*sin(theta)*sin(theta)*cos(phi)*sin(phi)','S0*sin(theta)*cos(phi)*cos(theta)','S0*(sin(theta)*sin(theta)*sin(phi)*sin(phi) - 0.5)', 'S0*sin(theta)*sin(phi)*cos(theta)'), theta = theta,phi = phi, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1, q_target2, q_target3, q_target4 = q_target_proj.split()

elif d == 2 and target_geometry == "circle_planar":
    theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
    q_target = Expression(('S0*(sin(theta)*sin(theta)-0.5)', '-S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1 = q_target_proj.split()

elif d == 2 and target_geometry == "defect":
    q_target = Expression(('0', '0'), degree = 1)
    for i, q in enumerate(config['top_charge']):
        position_defect = config['position_defect'][i]
        # create a subdomain of radius 0.5 around the defect with label i
        tol = 1E-12
        class Charged(SubDomain):
            def inside(self, x, on_boundary):
                r = ((x[0]-position_defect[0])**2 + (x[1]-position_defect[1])**2)**0.5
                return r <= 0.45 + tol
        Charged_Domain = Charged()
        Charged_Domain.mark(domains, i)
        theta = Expression('q*atan2((x[1]-y0),(x[0]-x0))', degree = 1, q = q, x0 = position_defect[0], y0 = position_defect[1])
        q_target += Expression(('sqrt(pow(x[0]-x0,2) + (pow(x[1]-y0,2))) < 0.5 ? S0*(cos(theta)*cos(theta)-0.5) : 0', 'sqrt(pow(x[0]-x0,2) + (pow(x[1]-y0,2))) < 0.5 ? S0*sin(theta)*cos(theta):0'), theta = theta, S0 = S0, x0 = position_defect[0], y0 = position_defect[1], degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1 = q_target_proj.split()

elif d == 3 and target_geometry == "pseudoChiral":
    q_target = Expression(('S0*(cos(phi)*cos(phi)-1/3)', 'S0*sin(phi)*cos(phi)', '0', 'S0*sin(phi)*sin(phi) - 1/3', '0'), phi = Expression('pi/4*0.2*x[2]', degree = 1) , S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1, q_target2, q_target3, q_target4 = q_target_proj.split()
else:
    q_target = Expression(('0', '0'), degree = 1)
    if d == 3:
        q_target = Expression(('0', '0', '0', '0', '0'), degree = 1)
    q_target_proj = project(q_target, W)
    if d == 2:
        q_target0, q_target1 = q_target_proj.split()
    else:
        q_target0, q_target1, q_target2, q_target3, q_target4 = q_target_proj.split() 
if d == 2:
    File(save_dir + '/target.pvd') << q_target_proj # Save the target Q
else: 
    xdmffile_target = XDMFFile(save_dir + '/target.xdmf')
    xdmffile_target.write(q_target_proj,0)


target_orientation, target_S = compute_orientation(q_target, mesh, d)
File(save_dir + '/target_director.pvd') << target_orientation # Save the target director field
File(save_dir + '/target_S.pvd') << target_S # Save the target scalar order parameter field
# Define the objective function

# The objective function is the sum of the squared difference between the state variable q_ and the target q_target,
# the volume constraint, the center of mass constraints 
# and a regularization term for the mesh quality.
# TODO Subtract forms
objective_main = ((q_[0]-q_target0)**2 + (q_[1]-q_target1)**2)*dx((0,1))  
if d == 3:
    objective_main += ((q_[2]-q_target2)**2 + (q_[3]-q_target3)**2 + (q_[4]-q_target4)**2)*dx
objective_volume = ((current_volume - Vol0)**2)/current_volume*dx
objective_center_of_mass = ((assemble(x*dx)-c_x0)**2 + (assemble(y*dx)-c_y0)**2)*dx
objective_boundary_length = (current_boundary_length - Boundary0)**2 /current_volume*dx
objective_mesh_quality = (1/(CellVolume(mesh)+Constant(eps_meshquality))**2)/current_volume*dx
objective_max_cell_edge_length = (MaxCellEdgeLength(mesh)**2)/current_volume*dx

objective = objective_main + k_meshquality*objective_mesh_quality + k_vol*objective_volume + k_com*objective_center_of_mass + k_boundarylength*objective_boundary_length + k_edge*objective_max_cell_edge_length

Lagrangian = objective + Dv_Energy
forwardPDE = derivative(Lagrangian, p_, TestFunction(W)) 
adjointPDE = replace(derivative(Lagrangian, q_, TestFunction(W)), {p_: TrialFunction(W)})
forwardJacobian = derivative(forwardPDE, q_)
adjointJacobian = derivative(adjointPDE, p_)


initial_guess = compute_initial_guess(mesh, S0, boundaries, surf_markers, finite_element, finite_element_degree, d, config['anchoring'])
assign(q_, initial_guess)
assign(p_, initial_guess)
if d == 2:
    File(save_dir + '/initialguess.pvd') << q_
else: 
    xdmffile_initialguess = XDMFFile(save_dir + '/initialguess.xdmf')
    xdmffile_initialguess.write(q_,0)
initial_orientation, initial_S = compute_orientation(q_, mesh, d)
File(save_dir + '/initial_director.pvd') << initial_orientation # Save the initial director field
File(save_dir + '/initial_S.pvd') << initial_S

# Define functions which are going to hold the solution of the shape gradient,
# as well as a scaled version of the negative shape gradient.
shapeGradient = Function(S)
displacement = Function(S)

# Setup the elasticity inner product.
delta = Constant(config['delta']) # Regularization parameter for the elasticity inner product
delta_beltrami = Constant(config['delta_beltrami']) # Regularization parameter for the elasticity inner product
if config['inner_product'] == "H1":
    displacementInnerProduct = inner(grad(TrialFunction(S)), grad(TestFunction(S))) * dx  + delta * inner(TrialFunction(S),TestFunction(S)) * dx
elif config['inner_product'] == "elasticity":
    displacementInnerProduct = inner(C(strain(TrialFunction(S))),strain(TestFunction(S))) * dx + delta * inner(TrialFunction(S),TestFunction(S)) * dx
else: 
    raise ValueError("Inner product not supported")

if config['tangential_smoothing']:
    displacementInnerProduct += delta_beltrami*inner(tang_grad(TrialFunction(S), normals), tang_grad(TestFunction(S), normals)) * ds_controlvariable

objective_values, shape_gradient_norms, rel_changes, objectives_main, objectives_meshquality, volumes, variances_radius = [], [], [], [], [], [], []
# # Plot the subdomains
# import matplotlib.pyplot as plt
# plt.figure()
# plot(domains, title="Subdomains")
# plt.savefig(save_dir + "/subdomains.png")
# plt.close()

# asd = assemble(1*dx((0,1)))
# print(asd)

# Run a shape gradient loop.
iteration = 0
alpha = alphaInit
if d == 3:
    xdmffile_results = XDMFFile(save_dir + '/results.xdmf')
    xdmffile_adjoints = XDMFFile(save_dir + '/adjoints.xdmf')

while iteration < maxIter:
    X = SpatialCoordinate(mesh)
    # Compute the current volume and boundary length of the mesh.
    current_volume.assign(assemble(1.0*dx))
    current_boundary_length.assign(assemble(1.0*ds_controlvariable))
    # Solve the forward PDE.
    # compute initial guess
    initial_guess = compute_initial_guess(mesh, S0, boundaries, surf_markers, finite_element, finite_element_degree, d, config['anchoring'])
    assign(q_, initial_guess)
    assign(p_, initial_guess)

    solveMultRelaxation([[0.3,1e-3],[0.4,1.0e-5], [1.0,1e-8]], forwardPDE,0, q_, None, forwardJacobian, ffc_options)
    
    J = assemble(objective)

    # Solve the adjoint PDE.
    solve(lhs(adjointPDE) == rhs(adjointPDE), p_, None)

    # Evaluate the shape derivative. 
    nh = compute_normals(mesh)
    dJds_main = derivative(objective + Dv_Energy, X, TestFunction(S))  
    dJdS_com = derivative(objective_center_of_mass, X, TestFunction(S))
    dJdS_meshquality = (div(TestFunction(S)) / (CellVolume(mesh) + Constant(eps_meshquality))) / current_volume * dx

    #  Warning: Check for correct implementation!!!
    dJds_volume = 2*(current_volume - Vol0)*inner(TestFunction(S), normals)/current_volume*ds_controlvariable
    dJds_boundary_length = 2*(current_boundary_length - Boundary0)*inner(TestFunction(S),normals)*div(nh)/current_volume*ds_controlvariable
    dJdS_max_cell_edge_length = (div(TestFunction(S)) * MaxCellEdgeLength(mesh))/current_volume * dx

    dJds = dJds_main

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

    # Use only the normal component of this boundary vector as BC
    N_vec = CGNormal(mesh)
    dJds_1_n = project(inner(dJds_1, N_vec) * N_vec, S)
    File(save_dir + f"/dJds_boundary_{iteration}.pvd") << dJds_1
    File(save_dir + f"/dJds_1_normal_{iteration}.pvd") << dJds_1_n
    File(save_dir + f"/N_vec_{iteration}.pvd") << N_vec
    # bc_outer = DirichletBC(S, dJds_1_n, 'on_boundary')

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
    
    M, S_param = compute_orientation(q_, mesh, d)
    File(save_dir + '/directorfield-{0:03d}.pvd'.format(iteration)) << M
    File(save_dir + '/scalarorderparameter-{0:03d}.pvd'.format(iteration)) << S_param
    neg_shape_gradient = Function(S)
    neg_shape_gradient.assign(-shapeGradient)
    File(save_dir + f"/neg_shape_gradient_{iteration}.pvd") << neg_shape_gradient

    # Printing information
    terms = [
    ("objective_main", objective_main, dJds_main, 1.0),
    ("objective_com", objective_center_of_mass, dJdS_com, float(k_com.values()[0])),
    # Add more terms as needed
        ]
    print(f'Iteration {iteration}')
    print_objective_and_gradient_terms(terms)
    objective_values.append(J)
    shape_gradient_norms.append(normShapeGradient2)
    objectives_main.append(assemble(objective_main))
    objectives_meshquality.append(assemble(objective_mesh_quality))
    # variances_radius.append(variance_radius(mesh, surf_markers, boundaries, dx))
    volumes.append(current_volume.values()[0])
    # Store the mesh associated with the current iterate, as well as its objective value.
    referenceMeshCoordinates = mesh.coordinates().copy()

    # Begin Armijo line search.
    lineSearchSuccessful = False
    sub_iteration = 0
    # if iteration == 0:
    #     alpha = 1/sqrt(normShapeGradient2)
    # else:
    #     alpha = alpha*shape_gradient_norms[-2]/shape_gradient_norms[-1]
    
    while (lineSearchSuccessful == False) and (alpha > alphaMin):
        # Assign the mesh displacement vector field.
        mesh.coordinates()[:] = referenceMeshCoordinates
        displacement.assign(- alpha * shapeGradient)

        # Update the mesh by adding the displacement to the reference mesh.
        ALE.move(mesh, displacement)
        File(save_dir + f"/mesh_iter{iteration}_sub{sub_iteration}.pvd") << mesh
        # Update the current volume and boundary length of the mesh.
        current_volume.assign(assemble(1.0*dx))
        current_boundary_length.assign(assemble(1.0*ds_controlvariable))

        # Solve the forward PDE.
        assign(q_, compute_initial_guess(mesh, S0, boundaries, surf_markers, finite_element, finite_element_degree, d, config['anchoring']))
        # Solve the forward PDE with the updated mesh.
        solveMultRelaxation([[0.3,1e-3],[0.4,1.0e-5], [1.0,1e-8]], forwardPDE,0, q_, None, forwardJacobian, ffc_options)

        trialJ = assemble(objective)


        # Evaluate the Armijo condition and reduce the step size if necessary.
        if (trialJ <= J - sigma * alpha * normShapeGradient2):
            lineSearchSuccessful = True
        else:
            alpha = beta * alpha

        # Write debugging information to a file
        with open(save_dir + '/debugging.txt', 'a') as debug_file:
            debug_file.write(f"It.: {iteration}\t")
            debug_file.write(f"Alpha (step size): {alpha:9.2e}\t")
            debug_file.write(f"Obj. (J): {J:9.2e}\t")
            debug_file.write(f"Trial Obj. (trialJ): {trialJ:9.2e}\t")
            debug_file.write(f"Norm of Shape Gradient squared: {normShapeGradient2:12.2e}\t")
            debug_file.write(f"Armijo Condition: {lineSearchSuccessful}\n")
            debug_file.write("-" * 40 + "\n")

        # Increment the sub-iteration counter.
        sub_iteration += 1
    # Occasionally display some information.

    # Reset the step size for the next iteration.
    alpha = min(alphaInit, 1000* alpha / beta)
    # Set the initial step size.
    if iteration > 0:
        rel_change = abs(objective_values[-1] - objective_values[-2]) / (abs(objective_values[-2]) + 1e-12)
        rel_changes.append(rel_change)
        if rel_change < float(config['rel_change_stopping_value']):
            print(f"Stopping: Relative change in objective ({rel_change:.2e}) is below threshold.")
            break

    # Increment the iteration counter.
    iteration += 1

write_objective_terms_to_file(save_dir, {'objective_values': objective_values, 'objectives_main': objectives_main, 'objectives_meshquality': objectives_meshquality, 'volumes': volumes, 'variance radii': variances_radius})
plotResults(save_dir, objective_values, shape_gradient_norms, rel_changes, objectives_main, objectives_meshquality, k_meshquality)




# TODO Symmetric Mesh
# Generate quarter mesh
# or only computer on quarter mesh and mirror the result

# Penalty length, Gauß divergence

# TODO Nabla grad of the normals in DG 0, gives 0?

# TODO Hand derive volume constraint maybe?

# TODO 1/2 defects
#  Cholersteric defect

# Novel W1,infty inner product
# Paper: https://arxiv.org/pdf/2103.13857 

# Stephan Schmidt:

# Solve the boundary smoothing problem use it as a DirichletBC with no source term for the volume defomration 

# Paper: Blind of the physics: Schmidt Schütte, Walther

# Mitoush Version can derive through

# Topology Derivative
# Paper: https://arxiv.org/abs/2303.15070
# Paper: https://arxiv.org/pdf/2411.19421
# Paper: https://arxiv.org/abs/2307.12444



# Projective gradient descent

# For DG0, inculde Jump terms? Over the skeleton dS int over jump(u) jump(u)

# Lowlevel test, MWE
# Inhomogenous mesh, finer mesh in center
# Adaptive scheme, a posteriori refine mesh

# Custom evaluation routine; Userdefined evaluation (search online), slow. Start with interpolation (maybe to higher (second) order)?