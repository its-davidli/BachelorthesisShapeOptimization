from fenics import *
from dolfin import *
from dolfin_adjoint import *
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
os.system("mkdir "+" "+save_dir + "/Figures_and_Data")

# Save the configuration file in the save directory
os.system("cp config.yaml " + save_dir + "/Figures_and_Data/config.yaml")
os.system("cp "+ __file__ + " " + save_dir + "/Figures_and_Data/code.py")
os.system("cp "+ 'DeliverableShapeOptimizationLCsMethods.py' + " " + save_dir + "/Figures_and_Data/methods.py")

# Folder for initial mesh files
subfolder = config['subfolder']
geom_folder = config['geom_folder']
mesh_name = config['mesh_name']

# Clear debugging file
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

# Visualize boundary subdomain with subdomain_ids=surf_markers, for visual verification of anchoring boundary
# 0 corresponds to non-anchoring boundary facets, 1 to anchoring boundary facets
boundary_marked = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_marked.set_all(0)
for marker in surf_markers_anchoring:
    for facet in facets(mesh):
        if boundaries[facet.index()] == marker:
            boundary_marked[facet.index()] = 1
File(save_dir + "/boundary_subdomain_anchoring.pvd") << boundary_marked

# Visualize boundary subdomain with subdomain_ids=surf_markers_moving, for visual verification of moving boundary
# 0 corresponds to non-moving boundary facets, 1 to moving boundary facets
boundary_marked = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_marked.set_all(0)
for marker in surf_markers_moving:
    for facet in facets(mesh):
        if boundaries[facet.index()] == marker:
            boundary_marked[facet.index()] = 1
File(save_dir + "/boundary_subdomain_moving.pvd") << boundary_marked

# Read physical parameters from config
n_indep_comp = int(d*(d+1)/2 - 1) #Number of independent components for a symmetric traceless tensor)
L_c = float(config["L_c"])
a_B = float(config["a_B"])
if d==2: S0 = np.sqrt(2.0*a_B) # Corresponding LdG S eigenvalue
if d==3: S0 = 1/4*(math.sqrt(24.0*a_B + 1)+1) # Corresponding LdG S eigenvalue

# Set default algorithmic parameters (Armijo line search)
maxIter = int(config['maxIter'])        # Maximum number of shape gradient steps.
sigma = float(config['sigma'])          # Armijo line search slope parameter
beta = float(config['beta'])            # Armijo backtracking parameter
alphaInit = float(config['alphaInit'])  # Initial step size for Armijo line search
alphaMin = float(config['alphaMin'])    # Fallback minimal step size to terminate Armijo line search.
tol = 1E-12                             # General tolerance for floating point comparisons

k_bc = Constant(float(config['k_bc']))  # Penalty constant for the (weak) anchoring term in the LdG energy functional 

# Define data and auxiliary functions for the elasticity inner product.
if config['inner_product'] == "elasticity":
    E = Constant(float(config['E'])) # Young's modulus
    nu = Constant(float(config['nu'])) # Poisson's ratio
    lmbda = nu*E/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))
    if d==2:
        lmbda = 2*mu*lmbda/(lmbda + 2*mu) # Plane stress condition
    def strain(u): return sym(nabla_grad(u))
    # def strain(u): return sym(nabla_grad(u)) - (1.0/3.0)*tr(sym(nabla_grad(u)))*Identity(d)
    def C(epsilon): return 2* mu * epsilon + lmbda * tr(epsilon) * Identity(d)

x = Expression('x[0]', degree = 1)
y = Expression('x[1]', degree = 1)
if d==3: z = Expression('x[2]', degree = 1)

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
elif config['anchoring'] == 'planar':
    if d == 2:
        tangentials = as_vector([-normals[1], normals[0]])
        Q_b = (S0/2.0)*(d*outer(tangentials, tangentials)-Identity(d)) 
    else:
        raise NotImplementedError("Tangential vector computation is only implemented for 2D meshes.")


Energy = (LG_energy(q_,a_B,L_c, d))*dx + LG_boundaryPenalty(q_,Q_b,k_bc, d) * ds_anchoring
Dv_Energy= derivative(Energy, q_, p_)

# Define the target state variable q_target
target_geometry = config['target_geometry']
subdomainlist = [] # List of subdomain ids where the objective function is evaluated
if d == 2 and target_geometry == "circle":
    X = SpatialCoordinate(mesh)
    theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
    q_target = Expression(('S0*(cos(theta)*cos(theta)-0.5)', 'S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
elif d == 2 and target_geometry == "half_of_circle":
    X = SpatialCoordinate(mesh)
    theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
    q_target = Expression(('S0*(cos(theta)*cos(theta)-0.5)', 'S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    class Charged(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= 0
    Charged_Domain = Charged()
    Charged_Domain.mark(domains, 1)
    subdomainlist.append(1)

elif d ==2 and target_geometry == "interior_of_circle":
    X = SpatialCoordinate(mesh)
    theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
    q_target = Expression(('S0*(cos(theta)*cos(theta)-0.5)', 'S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    class Charged(SubDomain):
        def inside(self, x, on_boundary):
            return x[0]**2 + x[1]**2 <= 0.5
    Charged_Domain = Charged()
    Charged_Domain.mark(domains, 1)
    subdomainlist.append(1)

elif d == 2 and target_geometry == "uniform_vertical":
    X = SpatialCoordinate(mesh)
    q_target = Expression(('-S0*(0.5)', 'eps'), S0 = S0, eps=tol, degree = 1)
    q_target_proj = project(q_target, W)
    # class Charged(SubDomain):
    #     def inside(self, x, on_boundary):
    #         return x[1] <= 1.5 and x[1] >= -1.5
    # Charged_Domain = Charged()
    # Charged_Domain.mark(domains, 100)
    # subdomainlist.append(100)    
    q_target_proj = project(q_target, W)

elif d == 3 and target_geometry == "sphere":
    phi = Expression('atan2((x[1]),(x[0]))', degree = 1)
    theta = Expression('acos(x[2]/std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+1e-20))', degree = 1)
    q_target = Expression(('S0*(cos(phi)*cos(phi)*sin(theta)*sin(theta)-0.5)','S0*sin(theta)*sin(theta)*cos(phi)*sin(phi)','S0*sin(theta)*cos(phi)*cos(theta)','S0*(sin(theta)*sin(theta)*sin(phi)*sin(phi) - 0.5)', 'S0*sin(theta)*sin(phi)*cos(theta)'), theta = theta,phi = phi, S0 = S0, degree = 1)
    q_target_proj = project(q_target, W)
    # q_target0, q_target1, q_target2, q_target3, q_target4 = q_target_proj.split()

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
        class Charged(SubDomain):
            def inside(self, x, on_boundary):
                r = ((x[0]-position_defect[0])**2 + (x[1]-position_defect[1])**2)**0.5
                return r <= 0.45 + tol
        Charged_Domain = Charged()
        Charged_Domain.mark(domains, i)
        subdomainlist.append(i)
        theta = Expression('q*(atan2((x[1]-y0),(x[0]-x0)+ (-ind + 1)* pi))', degree = 1, ind = i, q = q, x0 = position_defect[0], y0 = position_defect[1])
        q_target += Expression(('sqrt(pow(x[0]-x0,2) + (pow(x[1]-y0,2))) < 0.5 ? S0*(cos(theta)*cos(theta)-0.5) : 0', 'sqrt(pow(x[0]-x0,2) + (pow(x[1]-y0,2))) < 0.5 ? S0*sin(theta)*cos(theta):0'), theta = theta, S0 = S0, x0 = position_defect[0], y0 = position_defect[1], degree = 1)
    q_target_proj = project(q_target, W)
    q_target0, q_target1 = q_target_proj.split()

elif d == 3 and target_geometry == "pseudoChiral":
    X = SpatialCoordinate(mesh)
    phi = pi/4*0.2*X[2]
    eps = Constant(1e-10)
    # q_target = as_vector((S0*(cos(phi)*cos(phi)-1/3), S0*sin(phi)*cos(phi), eps, S0*sin(phi)*sin(phi) - 1/3, eps))
    q_target = Expression(('S0*(cos(phi)*cos(phi)-1/3)', 'S0*sin(phi)*cos(phi)', 'eps', 'S0*sin(phi)*sin(phi) - 1/3', 'eps'), phi = Expression('pi/4*0.2*x[2]', degree = 1) , S0 = S0, eps= Expression('pow(10,-10)', degree = 1), degree = 1)
    q_target_proj = project(q_target, W)
    # q_target0, q_target1, q_target2, q_target3, q_target4 = q_target_proj.split()

elif d == 3 and target_geometry == "uniform_horizontal":
    X = SpatialCoordinate(mesh)
    azimuth_angle = Constant(np.pi/6)
    q_target = Expression(('S0*(cos(phi)*cos(phi)-1/3)', 'S0*sin(phi)*cos(phi)', 'eps', 'S0*(sin(phi)*sin(phi) - 1/3)', 'eps'), phi = azimuth_angle , S0 = S0, eps= tol, degree = 1)
    q_target_proj = project(q_target, W)

elif d == 3 and target_geometry == "saturnring_defect":
    X = SpatialCoordinate(mesh)
    radius_circle = 0.05
    phi = Expression('atan2((x[1]),(x[0]))', degree = 1)
    theta = Expression('acos(x[2]/std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+1e-20))', degree = 1)
    Q_h = Expression(('S0*(cos(phi)*cos(phi)*sin(theta)*sin(theta)-0.5)','S0*sin(theta)*sin(theta)*cos(phi)*sin(phi)','S0*sin(theta)*cos(phi)*cos(theta)','S0*(sin(theta)*sin(theta)*sin(phi)*sin(phi) - 0.5)', 'S0*sin(theta)*sin(phi)*cos(theta)'), theta = theta,phi = phi, S0 = S0, degree = 1)
    Q_inf = Expression(('-S0*(1.0/3.0)', 'eps', 'eps', '-S0*(1.0/3.0)', 'eps'), S0 = S0, eps= 0.0, degree = 1)
    rescaled_radius = Expression('sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) / R', R = radius_circle, degree = 1)
    w = Constant(radius_circle*k_bc/L_c)
    q_target  = w/(3+w)/rescaled_radius**3*Q_h + (1 - w/(1+w)/rescaled_radius)*Q_inf
    q_target_proj = project(q_target, W)


elif d == 2 and target_geometry == "saturnring_defect":
    radius_circle = 0.025
    phi = Expression('atan2((x[1]),(x[0]))', degree = 1)
    Q_h = Expression(('S0*(cos(phi)*cos(phi)-0.5)', 'S0*sin(phi)*cos(phi)'), phi = phi, S0 = S0, degree = 1)
    Q_inf = Expression(('-S0*(0.5)', 'eps'), S0 = S0, eps= 0.0, degree = 1)
    rescaled_radius = Expression('sqrt(x[0]*x[0] + x[1]*x[1]) / R', R = radius_circle, degree = 1)
    w = Constant(radius_circle*k_bc/L_c)
    q_target  = w/(3+w)/rescaled_radius**2*Q_h + (1 - w/(1+w)/rescaled_radius)*Q_inf
    class Charged(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] <= 0.0 - 1.5 *radius_circle
    Charged_Domain = Charged()
    Charged_Domain.mark(domains, 10)
    subdomainlist.append(10)
    q_target_proj = project(q_target, W)

else:
    raise ValueError("Target geometry not supported.")

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
if not subdomainlist:
    objective_main = (dot(q_ - q_target, q_ - q_target))/(assemble(1*dx))*dx
else:
    objective_main = (dot(q_ - q_target, q_ - q_target))/(assemble(1*dx(subdomainlist[0])))*dx(subdomainlist[0])
    # Visualize subdomain where objective is evaluated
    area_marked = MeshFunction("size_t", mesh, mesh.topology().dim())
    area_marked.set_all(0)
    for facet in facets(mesh):
        if domains[facet.index()] == subdomainlist[0]:
            area_marked[facet.index()] = 1
    File(save_dir + "/area_obj_eval.pvd") << area_marked


# Combine objective with the LdG energy to form the Lagrangian of the shape optimization problem and set up adjoint and forward problems
objective = objective_main
Lagrangian = objective + Dv_Energy
forwardPDE = derivative(Lagrangian, p_, TestFunction(W)) 
adjointPDE = replace(derivative(Lagrangian, q_, TestFunction(W)), {p_: TrialFunction(W)})
forwardJacobian = derivative(forwardPDE, q_)
adjointJacobian = derivative(adjointPDE, p_)


# Compute and save initial guess for the state variable q
initial_guess = compute_initial_guess(mesh, config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring'])
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
if config['inner_product'] == "H1":
    displacementInnerProduct = inner(grad(TrialFunction(S)), grad(TestFunction(S))) * dx  
elif config['inner_product'] == "elasticity":
    displacementInnerProduct = inner(C(strain(TrialFunction(S))),strain(TestFunction(S))) * dx 
else: 
    raise ValueError("Inner product not supported")

delta = Constant(config['delta']) # Regularization parameter for the elasticity inner product
delta_beltrami = Constant(config['delta_beltrami']) # Regularization parameter for the elasticity inner product
displacementInnerProduct += delta*inner(TrialFunction(S), TestFunction(S)) * dx
if config['tangential_smoothing']:
    displacementInnerProduct += delta_beltrami*inner(tang_grad(TrialFunction(S), normals), tang_grad(TestFunction(S), normals)) * ds_controlvariable

objective_values, alphas, shape_gradient_norms, rel_changes, objectives_main, objectives_meshquality, volumes, variances_radius, radii, center_of_masses = [], [], [], [], [], [], [], [], [], []

# Define boundary conditions for the shape gradient, if needed
bc_shapegradient = []
if config['boundary_conditions_shapegradient'] == 'None':
    pass

elif d == 3 and config['boundary_conditions_shapegradient'] == 'fixed_bottom':
    def boundary(x):
        return near(x[2], 0, tol)
    bc_shapegradient += [DirichletBC(S, Expression(('0','0','0'), degree = 1), boundary)]

elif d==2 and config['boundary_conditions_shapegradient'] == 'fixed_sides':
    bc_shapegradient += [DirichletBC(S, Expression(('0','0'), degree = 1), boundaries, 2)]


elif d==2 and config['boundary_conditions_shapegradient'] == 'fixed_square':
    def boundary(x):
        return near(abs(x[0]), 0.5, tol) or near(abs(x[1]), 0.5, tol)
    bc_shapegradient += [DirichletBC(S, Expression(('0','0'), degree = 1), boundary)]

elif d==3 and config['boundary_conditions_shapegradient'] == 'fixed_sides':
    bc_shapegradient += [DirichletBC(S, Expression(('0','0','0'), degree = 1), boundaries,config['boundary_conditions_shapegradient_markers'][0])]

else:
    raise ValueError("Boundary condition for shape gradient not supported")

# Set up XDMF files for saving results in 3D
if d == 3:
    xdmffile_results = XDMFFile(save_dir + '/results.xdmf')
    xdmffile_adjoints = XDMFFile(save_dir + '/adjoints.xdmf')

# Run a shape gradient loop.
iteration = 0
while iteration < maxIter:
    X = SpatialCoordinate(mesh)

    # Solve the forward PDE.
    # compute initial guess
    initial_guess =  compute_initial_guess(mesh, config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring'])
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

    # Use only the normal component of this boundary vector as BC
    N_vec = CGNormal(mesh)
    dJds_1_n = project(inner(dJds_1, N_vec) * N_vec, S)
    File(save_dir + f"/dJds_rep_{iteration}.pvd") << dJds_rep
    File(save_dir + f"/dJds_boundary_{iteration}.pvd") << dJds_1
    File(save_dir + f"/dJds_1_normal_{iteration}.pvd") << dJds_1_n
    File(save_dir + f"/N_vec_{iteration}.pvd") << N_vec

    if config['use_normal_dJds_projection']:
        solve(displacementInnerProduct == inner(dJds_1_n,TestFunction(S))*dx, shapeGradient, bc_shapegradient)
    else:
        solve(displacementInnerProduct == dJds, shapeGradient, bc_shapegradient)


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
    # Add more terms as needed
        ]
    print(f'Iteration {iteration}')
    print_objective_and_gradient_terms(terms)
    objective_values.append(J)
    shape_gradient_norms.append(normShapeGradient2)
    objectives_main.append(assemble(objective_main))
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
        if config['stepsize_method'] == "herzog_meshquality":
            alpha = 1/sqrt(normShapeGradient2)
        if config['stepsize_method'] == "constant":
            alpha = alphaInit
    else:
        if config['stepsize_method'] == "herzog_meshquality":
            alpha = alpha*sqrt(shape_gradient_norms[-2]/normShapeGradient2)
            if alpha * sqrt(normShapeGradient2) < 1e-3:
                alpha = 1/sqrt(normShapeGradient2)
        if config['stepsize_method'] == "constant":
            alpha = min(alphaInit, 1000* alpha / beta)
    # alpha = 1.25 * alpha
    while (lineSearchSuccessful == False) and (alpha > alphaMin):
        # Assign the mesh displacement vector field.
        mesh.coordinates()[:] = referenceMeshCoordinates
        displacement.assign(- alpha * shapeGradient)

        # Update the mesh by adding the displacement to the reference mesh.
        ALE.move(mesh, displacement)
        File(save_dir + f"/mesh_iter{iteration}_sub{sub_iteration}.pvd") << mesh
        
        # Solve the forward PDE.
        assign(q_, compute_initial_guess(mesh, config['initial_guess'], S0, boundaries, surf_markers_anchoring, finite_element, finite_element_degree, d, config['anchoring']))
        # Solve the forward PDE with the updated mesh.
        solveMultRelaxation(config['forward_solver_relaxations'], forwardPDE,0, q_, None, forwardJacobian, ffc_options)

        trialJ = assemble(objective)


        # Evaluate the Armijo condition and reduce the step size if necessary.
        if (trialJ <= J - sigma * alpha * normShapeGradient2):
            lineSearchSuccessful = True
            alphas.append(alpha)
        
        # Write debugging information to a file
        with open(save_dir + '/Figures_and_Data/debugging.txt', 'a') as debug_file:
            debug_file.write(f"It.: {iteration}\t")
            debug_file.write(f"Alpha (step size): {alpha:9.2e}\t")
            debug_file.write(f"Obj. (J): {J:9.2e}\t")
            debug_file.write(f"Trial Obj. (trialJ): {trialJ:9.2e}\t")
            debug_file.write(f"Norm of Shape Gradient squared: {normShapeGradient2:12.2e}\t")
            debug_file.write(f"Armijo Condition: {lineSearchSuccessful}\n")
            debug_file.write("-" * 40 + "\n")

        alpha = beta * alpha

        # Increment the sub-iteration counter.
        sub_iteration += 1
    # Occasionally display some information.


    if iteration > 0:
        rel_change = abs(objective_values[-1] - objective_values[-2]) / (abs(objective_values[-2]) + 1e-12)
        rel_changes.append(rel_change)
        if rel_change < float(config['rel_change_stopping_value']):
            print(f"Stopping: Relative change in objective ({rel_change:.2e}) is below threshold.")
            break
        if trialJ> objective_values[-1]:
            print(f"Stopping: Objective functional increased: {objective_values[-1]} > {objective_values[-2]}.")
            break

    # Increment the iteration counter.
    write_objective_terms_to_file(save_dir, {'objective_values': objective_values, 'shape_gradient_norms_squared': shape_gradient_norms, 'alphas': alphas})
    plotResults(save_dir, objective_values, shape_gradient_norms, rel_changes)
    plotGeometricalInformation(save_dir, radii,variances_radius, center_of_masses)
    iteration += 1