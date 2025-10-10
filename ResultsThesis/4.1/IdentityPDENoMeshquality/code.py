from fenics import *
from dolfin import *
from dolfin_adjoint import *
import numpy as np
from ufl_legacy import nabla_div, nabla_grad, VectorElement, FiniteElement, MixedElement, split, atan_2, replace, Jacobian, log
import sys
import os
from IPython import embed
import matplotlib.pyplot as plt
import yaml # For reading configuration files
set_log_level(LogLevel.WARNING)
from Methods2DInnerMesh import (
    Stick,
    Ellipse,
    Circle,
    PiecewiseSource,
    compute_and_save_electric_field,
    print_objective_and_gradient_terms,
    mark_interface_facets,
    CGNormal,
)
 
# Load configuration from YAML file
with open("config.yaml", 'r') as stream:
    config = yaml.safe_load(stream)

#Create folder to save results
save_fold = config['save_fold']
save_subdir = sys.argv[1] 
save_dir = "./"+save_fold + "/" + save_subdir
os.system("mkdir "+" "+save_dir)
os.system("mkdir "+" "+save_dir + "/Figures")

# Save the configuration file in the save directory
os.system("cp config.yaml " + save_dir + "/config.yaml")
os.system("cp "+ __file__ + " " + save_dir + "/code.py")

# Folder for initial mesh files
subfolder = config['subfolder']
geom_folder = config['geom_folder']
mesh_name = config['mesh_name']

# Clear debugging file
open(save_dir + '/debugging.txt', 'w').close

finite_element = config['finite_element']
finite_element_degree = config['finite_element_degree']

surf_markers = config['surf_markers']

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
k_boundary = Constant(float(config['k_boundary']))                          # Penalty constant for the boundary length constraint
k_meshquality = Constant(float(config['k_meshquality']))                    # Penalty constant for the mesh quality penalty
eps_meshquality = Constant(float(config['eps_meshquality']))                # Regularize the meshquality term

# Define data and auxiliary functions for the elasticity inner product.
E = Constant(float(config['E']))            # Young's modulus
nu = Constant(float(config['nu']))          # Poisson's ratio
lmbda1 = nu*E/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))
lmbda = 2*lmbda1*mu/(lmbda1 + 2*mu)
def strain(u): return sym(nabla_grad(u))
def C(epsilon): return lmbda * epsilon + 2 * mu * tr(epsilon) * Identity(d)

#Import .xdmf geometry file into mesh + Mesh value collection for boundary conditions
mesh = Mesh()

mvc_subdomain = MeshValueCollection("size_t", mesh, mesh.topology().dim())
mvc_boundaries = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)

with XDMFFile(MPI.comm_world, '../'+geom_folder+subfolder+"/"+mesh_name+"_2D.xdmf") as xdmf_infile:
    xdmf_infile.read(mesh)
    xdmf_infile.read(mvc_subdomain, "")

with XDMFFile(MPI.comm_world, '../'+geom_folder+subfolder+"/"+mesh_name+"_1D.xdmf") as xdmf_infile:
    xdmf_infile.read(mvc_boundaries, "")

domains = cpp.mesh.MeshFunctionSizet(mesh, mvc_subdomain)
boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc_boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=domains)
q_degree = 5
dx = dx(metadata={'quadrature_degree': q_degree}) 
# ds = Measure('ds', domain=mesh, subdomain_data=boundaries, subdomain_id = 1)
dS = Measure('dS', domain=mesh, subdomain_data=boundaries, subdomain_id=2)
uncharged_marker = 1
charged_marker = 2

W = FunctionSpace(mesh, finite_element, finite_element_degree)
W0 = FunctionSpace(mesh,'DG',0 )


# define forward problem
q_ = Function(W)
p_ = Function(W)
d = q_.geometric_dimension() 
Q = 1
f = Function(W0)
# forwardPDEConstraint = inner(nabla_grad(q_), nabla_grad(p_)) * dx  - f * p_ * dx
# forwardPDEConstraint = inner(nabla_grad(q_), nabla_grad(p_)) * dx + inner(q_, p_) * dx - f * p_ * dx
forwardPDEConstraint = q_ *p_ * dx  - f * p_ * dx

# Define boundary conditions
# bc_outer = DirichletBC(W, Expression('std::log(1/sqrt(x[0]*x[0] + x[1]*x[1]))', degree=1), boundaries, uncharged_marker)
bc_outer = DirichletBC(W, Constant(0), boundaries, uncharged_marker)
bcs = [bc_outer]
# bcs = []

# Define the objective
target_potential = Function(W)
objective_main = (q_ - target_potential)**2 * dx
# objective_main = (-grad(q_) - compute_and_save_electric_field(target_potential))**2 * dx
objective_mesh_quality = (1/(CellVolume(mesh)+Constant(eps_meshquality))**2)*dx
objective = objective_main + k_meshquality * objective_mesh_quality
Lagrangian = objective + forwardPDEConstraint

forwardPDE = replace(derivative(Lagrangian, p_, TestFunction(W)), {q_: TrialFunction(W)})
adjointPDE = replace(derivative(Lagrangian, q_, TestFunction(W)), {p_: TrialFunction(W)})
forwardJacobian = derivative(forwardPDE, q_)
adjointJacobian = derivative(adjointPDE, p_)

#Calculate Target potential
f.assign(project(Ellipse(a=2, b=1.5, Q= Q, degree=0), W0))
File(save_dir + "/target_f.pvd") << f

solve(lhs(forwardPDE) == rhs(forwardPDE), q_, bcs = bcs, solver_parameters={"linear_solver": "mumps"})
target_potential.assign(q_)
File(save_dir + "/target_potential.pvd") << target_potential
target_field = compute_and_save_electric_field(target_potential, name = "target_electric_field.pvd", save_dir=save_dir, mesh=mesh)

# (Re-)set the source 
f.assign(project(PiecewiseSource(domains, Q =Q, charged_marker=charged_marker, degree=0), W0))

# Functionspaces for the deformation
S = VectorFunctionSpace(mesh, "CG", 2)
s = Function(S)

ALE.move(mesh, s)
# Define functions which are going to hold the solution of the shape gradient,
# as well as a scaled version of the negative shape gradient.
shapeGradient = Function(S)
displacement = Function(S)

normals = FacetNormal(mesh)
# Setup the (regularized) elasticity inner product
delta = 0.2
def tang_grad(u, normals=normals):
    return grad(u) - outer(dot(grad(u), normals), normals)

# displacementInnerProduct = inner(C(strain(TrialFunction(S))),strain(TestFunction(S))) * dx + delta * inner(TrialFunction(S),TestFunction(S)) * dx
displacementInnerProduct = inner(grad(TrialFunction(S)), grad(TestFunction(S))) * dx + delta * inner(TrialFunction(S), TestFunction(S)) * dx
# Apply zero displacement Dirichlet BC on the outer boundary to prevent mesh movement there
bc_shape_outer = DirichletBC(S, Constant((0.0,)*d), boundaries, uncharged_marker)
shape_bcs = [bc_shape_outer]
# Run a shape gradient loop.
iteration = 0
alpha = alphaInit
weight_change_counter = 0
objective_values = []
shape_gradient_norms = []
rel_changes = []
while iteration < maxIter:
    
    X = SpatialCoordinate(mesh)

    # Solve the forward problem
    solve(lhs(forwardPDE) == rhs(forwardPDE), q_, bcs = bcs, solver_parameters={"linear_solver": "mumps"})
    
    # Calculate the objective
    J = assemble(objective)

    # Solve the adjoint PDE
    solve(lhs(adjointPDE) == rhs(adjointPDE), p_, bcs=bcs, solver_parameters={"linear_solver": "mumps"})

    # Compute the shape gradient
    # Compute Shape Derivative
    dJds_main = derivative(objective_main + forwardPDEConstraint, X, TestFunction(S))
    dJdS_meshquality = (div(TestFunction(S)) / (CellVolume(mesh) + Constant(eps_meshquality))) * dx
    dJds = dJds_main + k_meshquality * dJdS_meshquality
    # --- Inner interface Riesz map with normal projection (adapted from v1 algorithm) ---
    # save facet markers
    # Riesz representation on the inner interface: find u in S s.t.
    # ∫_Γ <u, v> dS = dJds(v) for all v in S, with Γ the interface
    u_b = TrialFunction(S)
    v_b = TestFunction(S)
    a_b = inner(u_b('+'), v_b('+')) * dS
    L_b = dJds
    A_b = assemble(a_b, keep_diagonal=True)
    b_b = assemble(L_b)
    A_b.ident_zeros()
    dJds_1 = Function(S)
    solve(A_b, dJds_1.vector(), b_b)

    # Use only the normal component at the interface as RHS load (no volume projection)
    # Compute a CG1 normal field concentrated on the marked interface facets (tag 10)
    N_vec = CGNormal(mesh, FacetMarker=boundaries, Region=[2])
    File(save_dir + f"/dJds_interface_vector_{iteration}.pvd") << dJds_1
    L_shape = inner(inner(dJds_1('+'), N_vec('+')) * N_vec('+'), TestFunction(S)('+')) * dS
    File(save_dir + f"/CGNormal_{iteration}.pvd") << N_vec

    # Solve shape gradient equation with interface RHS, keep outer boundary fixed
    solve(displacementInnerProduct == dJds, shapeGradient, bcs=shape_bcs)
    # solve(displacementInnerProduct == L_shape, shapeGradient, bcs=shape_bcs)
    # Evaluate the squared norm of the shape gradient induced by the (regularized)
    # elasticity inner product.
    normShapeGradient2 = sum(shapeGradient.vector() * shapeGradient.vector())

    # Save the forward solution, electric field, source term, adjoint solution and negative shape gradient
    File(save_dir + f"/potential_{iteration}.pvd") << q_
    compute_and_save_electric_field(q_, name = f"electric_field_{iteration}.pvd", save_dir=save_dir, mesh=mesh)
    File(save_dir + f"/f_{iteration}.pvd") << f
    File(save_dir + f"/adjoint_solution_{iteration}.pvd") << p_
    neg_shape_gradient = Function(S)
    neg_shape_gradient.assign(-shapeGradient)
    File(save_dir + f"/neg_shape_gradient_{iteration}.pvd") << neg_shape_gradient

    # Printing information
    terms = [
    ("objective_main", objective_main, dJds_main, 1.0),
    ("objective_mesh_quality", objective_mesh_quality, dJdS_meshquality, float(k_meshquality.values()[0])),
    # Add more terms as needed
        ]
    print(f'Iteration {iteration}')
    print_objective_and_gradient_terms(terms)

    # Store the mesh associated with the current iterate, as well as its objective value.
    referenceMeshCoordinates = mesh.coordinates().copy()

    # Begin Armijo line search.
    lineSearchSuccessful = False
    sub_iteration = 0


    # Store initial values before line search
    objective_values.append(J)
    shape_gradient_norms.append(normShapeGradient2)
    while (not lineSearchSuccessful) and (alpha > alphaMin):
        # Assign the mesh displacement vector field.
        mesh.coordinates()[:] = referenceMeshCoordinates
        displacement.assign(-alpha * shapeGradient)
        # Update the mesh by adding the displacement to the reference mesh.
        ALE.move(mesh, displacement)    

        # # Save the source of the subiterations
        # if iteration == 2 :
        #     File(save_dir + f"/mesh_subiteration_{iteration}_{sub_iteration}.pvd") << f

        # Solve the forward PDE with the updated mesh.
        solve(lhs(forwardPDE) == rhs(forwardPDE), q_, bcs=bcs, solver_parameters={"linear_solver": "mumps"})

        trialJ = assemble(objective)
        # Evaluate the Armijo condition and reduce the step size if necessary.
        if (trialJ <= J - sigma * alpha * normShapeGradient2):
            lineSearchSuccessful = True
        else:
            alpha = beta * alpha
        trial_main_objective = assemble(objective_main)
        # Write debugging information to a file
        with open(save_dir + '/debugging.txt', 'a') as debug_file:
            debug_file.write(f"It.: {iteration}\t")
            debug_file.write(f"Alpha (step size): {alpha:9.2e}\t")
            debug_file.write(f"Obj. (J): {J:9.2e}\t")
            debug_file.write(f"Trial Obj. (trialJ): {trialJ:9.2e}\t")
            debug_file.write(f"Main Obj. (trial_main_objective): {trial_main_objective:9.2e}\t")
            debug_file.write(f"Norm of Shape Gradient squared: {normShapeGradient2:12.2e}\t")
            debug_file.write(f"Armijo Condition: {lineSearchSuccessful}\n")
            debug_file.write(f"Counter: {weight_change_counter}\n")
            debug_file.write("-" * 40 + "\n")

        sub_iteration += 1
    # Reduce mesh quality parameter if Armijo condition is not met
    if alpha < alphaMin*1e2 and weight_change_counter < 4:
        weight_change_counter += 1
        k_meshquality.assign(float(config['k_meshquality'])/(10**weight_change_counter))
    alpha = min(alphaInit, 1e5* alpha / beta)
        # Set the initial step size.
    # alpha = alphaInit
    # Stopping criterion: relative change in objective
    if iteration > 0:
        rel_change = abs(objective_values[-1] - objective_values[-2]) / (abs(objective_values[-2]) + 1e-12)
        rel_changes.append(rel_change)
        if rel_change < 1e-2:
            print(f"Stopping: Relative change in objective ({rel_change:.2e}) is below threshold.")
            break
    # Increment the iteration counter.
    iteration += 1

fig, ax1 = plt.subplots()

color1 = 'tab:blue'
ax1.set_xlabel("Iteration", fontsize=12)
ax1.set_ylabel("Objective Value", color=color1, fontsize=12)
ax1.plot(objective_values, label="Objective Value", color=color1)
# Plot with log scales on both axes
ax1.set_yscale('log')
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel("Shape Gradient Norm Squared", color=color2, fontsize=12)
ax2.plot(shape_gradient_norms, label="Shape Gradient Norm Squared", color=color2)
ax2.set_yscale('log')
ax2.tick_params(axis='y', labelcolor=color2)

fig.tight_layout()
plt.savefig(save_dir + "/Figures/objective_and_gradient_norms_log.png")
plt.close()

# Plot with linear scales on both axes
fig_lin, ax1_lin = plt.subplots()
color1_lin = 'tab:blue'
ax1_lin.set_xlabel("Iteration", fontsize=12)
ax1_lin.set_ylabel("Objective Value", color=color1_lin, fontsize=12)
ax1_lin.plot(objective_values, label="Objective Value", color=color1_lin)
ax1_lin.tick_params(axis='y', labelcolor=color1_lin)

ax2_lin = ax1_lin.twinx()
color2_lin = 'tab:red'
ax2_lin.set_ylabel("Shape Gradient Norm Squared", color=color2_lin, fontsize=12)
ax2_lin.plot(shape_gradient_norms, label="Shape Gradient Norm Squared", color=color2_lin)
ax2_lin.tick_params(axis='y', labelcolor=color2_lin)

fig_lin.tight_layout()
plt.savefig(save_dir + "/Figures/objective_and_gradient_norms_linear.png")
plt.close()

# plot rel_changes 
plt.figure()
plt.plot(rel_changes, label="Relative Change in Objective")
plt.xlabel("Iteration", fontsize=12)
plt.legend(fontsize=12)
plt.savefig(save_dir + "/Figures/rel_changes.png")
plt.close()

