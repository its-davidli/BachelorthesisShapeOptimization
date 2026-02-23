from fenics import *
from dolfin import *
from ufl_legacy import nabla_div, nabla_grad, VectorElement, FiniteElement, MixedElement, split, atan_2, replace, Jacobian,  tr, variable, shape, Max, sqrt, Min, det, Identity, max_value, min_value
import numpy as np
import math

def LG_energy(q,a_B,L_c, d): #Landau de Gennes energy functional with a single elastic constant L
    if d == 2: 
        Q = as_tensor(((q[0],q[1]),(q[1],-q[0])))
    elif d == 3: 
        Q = as_tensor([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], -q[0]-q[3]]])
    else: 
        raise ValueError("Dimension not supported")
    
    dQ = nabla_grad(Q)
    tQ2 = inner(Identity(d),dot(Q,Q))
    tQ3 = inner(Identity(d),dot(Q,dot(Q,Q)))
    ψ = - a_B*tQ2/2.0 - tQ3/3.0 + tQ2**2/4.0
    E = L_c*inner(dQ,dQ)/2.0 + ψ
    return E

def LG_boundaryPenalty(q, Q_b, k_bc, d): # Penalty Term in the LdG energy functional for the boundary condition
    if d == 2: 
        Q = as_tensor(((q[0],q[1]),(q[1],-q[0])))
    elif d == 3: 
        Q = as_tensor([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], -q[0]-q[3]]])
    else: 
        raise ValueError("Dimension not supported")
    return 0.5*k_bc*inner(Q - Q_b, Q - Q_b)

def compute_initial_guess(mesh, W, type, S0, boundaries, surf_markers, finite_element, finite_element_degree, d, anchoring): # Compute initial guess by minimizing the Frank energy functional

    if type == "uniform" and d == 2:
        return project(as_vector((-S0*(0.5), 1e-14)), W)

    if type == "uniform" and d == 3:
        # return project(as_vector((S0*(2/3),1e-14, 1e-14, -S0*(1/3), 1e-14)), W)
        phi = np.pi/6
        return project(as_vector((S0*(math.cos(phi)**2- 1/3),S0*math.cos(phi)*math.sin(phi), 1e-7, S0*(math.sin(phi)*math.sin(phi)-1/3), 1e-7)), W)
    if type == "uniform_in_z" and d == 3:
        return project(as_vector((-S0*(1/3),1e-14, 1e-14, -S0*(1/3), 1e-14)), W)

    if type == "Frank":
        # Define suitable vector, tensor and scalar spaces
        T = TensorElement(finite_element, mesh.ufl_cell(), finite_element_degree) # Needed for the projection of Q_Guess

        # Define the function space for the Frank Energy minimization
        V = VectorElement(finite_element, mesh.ufl_cell(), finite_element_degree)
        Vspace = FunctionSpace(mesh, V)

        # Compute the normal vector at the boundary for the minimization of the Frank energy
        normals = FacetNormal(mesh)
        if anchoring == 'planar':
            normals = as_vector([-normals[1], normals[0]])
        u = TrialFunction(Vspace)
        v = TestFunction(Vspace)
        a = inner(u, v) * ds
        l = inner(normals, v) * ds
        A = assemble(a, keep_diagonal=True)
        L = assemble(l)
        A.ident_zeros()
        n = Function(Vspace)
        solve(A, n.vector(), L)

        # Minimize Frank functional to use as initial guess in LdG minimization
        n_frank = TrialFunction(Vspace)
        n0 = Function(Vspace)
        ϕ = TestFunction(Vspace)
        zero_vector = Constant((0.0,)*d)
        a = inner(nabla_grad(n_frank), nabla_grad(ϕ)) * dx
        l = inner(zero_vector, ϕ) * dx
        bc_f = [DirichletBC(Vspace, n, boundaries, surf_marker) for surf_marker in surf_markers]
        solve(a == l, n0, bc_f)

        Q_guess = project(S0 * ((outer(n0, n0) / inner(n0, n0)) - Identity(d) / d), FunctionSpace(mesh, T))
        if d == 2:
            q_guess = [Q_guess.sub(0), Q_guess.sub(1)]
        elif d == 3:
            q_guess = [Q_guess.sub(0), Q_guess.sub(1), Q_guess.sub(2), Q_guess.sub(4), Q_guess.sub(5)]
        else:
            raise ValueError("Dimension not supported")
        return q_guess
    else:
        raise ValueError("Initial guess type not supported")

def solveMultRelaxation(parameters, lhs, rhs,solFunc, bcs, J, form_compiler_parameters): # Function to solve the PDE with multiple relaxation parameters
    for arr in parameters:
        solve(lhs == rhs,solFunc,bcs,J=J, form_compiler_parameters=form_compiler_parameters, solver_parameters={"newton_solver":
                                    {"relative_tolerance": arr[1],"absolute_tolerance": 1.0e-7, 'maximum_iterations': 500, "relaxation_parameter":arr[0], "linear_solver":"mumps"}})
        

def get_elasticity_operators(E, nu, d, inner_product):
    E = Constant(float(E))  # Young's modulus
    nu = Constant(float(nu))  # Poisson's ratio
    lmbda = nu * E / ((1 + nu) * (1 - 2 * nu)) # Lame's first parameter
    mu = E / (2 * (1 + nu)) # Lame's second parameter

    if d == 2:
        # Plane stress condition
        lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

    def strain(u):
        eps = sym(nabla_grad(u))
        if inner_product == 'elasticity_trace_free':
            return eps - (1.0 / d) * tr(eps) * Identity(d)
        else:
            return eps

    def C(epsilon):
        return 2 * mu * epsilon + lmbda * tr(epsilon) * Identity(d)

    return strain, C

def eigenstate3_legacy(A):
    # https://github.com/michalhabera/dolfiny/blob/master/src/dolfiny/invariants.py
    """Eigenvalues and eigenprojectors of the 3x3 (real-valued) tensor A.
    Provides the spectral decomposition A = sum_{a=0}^{2} λ_a * E_a
    with eigenvalues λ_a and their associated eigenprojectors E_a = n_a^R x n_a^L
    ordered by magnitude.
    The eigenprojectors of eigenvalues with multiplicity n are returned as 1/n-fold projector.

    Note: Tensor A must not have complex eigenvalues!
    """
    if shape(A) != (3, 3):
        raise RuntimeError(f"Tensor A of shape {shape(A)} != (3, 3) is not supported!")
    #
    eps = 1.0e-10
    #
    A = variable(A)
    #
    # --- determine eigenvalues λ0, λ1, λ2
    #
    # additively decompose: A = tr(A) / 3 * I + dev(A) = q * I + B
    q = tr(A) / 3
    B = A - q * Identity(3)
    # observe: det(λI - A) = 0  with shift  λ = q + ω --> det(ωI - B) = 0 = ω**3 - j * ω - b
    j = tr(B * B) / 2  # == -I2(B) for trace-free B, j < 0 indicates A has complex eigenvalues
    b = tr(B * B * B) / 3  # == I3(B) for trace-free B
    # solve: 0 = ω**3 - j * ω - b  by substitution  ω = p * cos(phi)
    #        0 = p**3 * cos**3(phi) - j * p * cos(phi) - b  | * 4 / p**3
    #        0 = 4 * cos**3(phi) - 3 * cos(phi) - 4 * b / p**3  | --> p := sqrt(j * 4 / 3)
    #        0 = cos(3 * phi) - 4 * b / p**3
    #        0 = cos(3 * phi) - r                  with  -1 <= r <= +1
    #    phi_k = [acos(r) + (k + 1) * 2 * pi] / 3  for  k = 0, 1, 2
    p = 2 / sqrt(3) * sqrt(j + eps**2)  # eps: MMM
    r = 4 * b / p**3
    r = max_value(min_value(r, +1 - eps), -1 + eps)  # eps: LMM, MMH
    phi = acos(r) / 3
    # sorted eigenvalues: λ0 <= λ1 <= λ2
    λ0 = q + p * cos(phi + 2 / 3 * pi)  # low
    λ1 = q + p * cos(phi + 4 / 3 * pi)  # middle
    λ2 = q + p * cos(phi)  # high
    #
    # --- determine eigenprojectors E0, E1, E2
    #
    E0 = diff(λ0, A).T
    E1 = diff(λ1, A).T
    E2 = diff(λ2, A).T



    return [λ0, λ1, λ2], [E0, E1, E2]

def compute_orientation(q, mesh, d):
    if d == 2:
        eps = 1e-10
        Q = as_tensor(((q[0],q[1]),(q[1],-q[0])))
        lmbda1 = sqrt(-1*det(Q))
        M = as_vector(( q[0] + lmbda1, 1* q[1]))/((q[0] + lmbda1)**2 + q[1]**2 + eps)** 0.5
        return project(M, VectorFunctionSpace(mesh, 'CG', 1)), project(2*lmbda1, FunctionSpace(mesh, 'DG', 0))
    elif d == 3:
        λ, E = eigenstate3_legacy(as_tensor([[q[0], q[1], q[2]],
                                                       [q[1], q[3], q[4]],
                                                       [q[2], q[4], -q[0]-q[3]]]))

    v = E[2] * as_vector([0, 0, 1]) 
    v = v / sqrt(inner(v, v))  # Normalize the eigenvector 
    return project(v, VectorFunctionSpace(mesh, 'CG', 1), form_compiler_parameters={'quadrature_degree': 2}),  project(3/2*λ[2], FunctionSpace(mesh, 'DG', 0),form_compiler_parameters={'quadrature_degree': 2})  # Return the principal eigenvector (v0 corresponds to the largest eigenvalue λ0)

def compute_normals(mesh):
    n = FacetNormal(mesh)
    V = VectorFunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(u,v)*ds
    l = inner(n, v)*ds
    A = assemble(a, keep_diagonal=True)
    L = assemble(l)

    A.ident_zeros()
    nh = Function(V)

    solve(A, nh.vector(), L)

    return nh

  #compute the CG1 outer normal of mesh in the Region labeled by FacetMarker

    #computes normal everywhere if nothing is specified

    #works for shells and full meshes

def CGNormal(mesh, FacetMarker = None, Region = None):


    MaxDim = mesh.geometry().dim()
    SpaceV = VectorFunctionSpace(mesh, "CG", 1, MaxDim)
    N = Function(SpaceV)
    N.rename("PointNormal", "")
    values = N.vector().get_local()
    v2d = vertex_to_dof_map(SpaceV)
    v2d = v2d.reshape((-1, MaxDim))
    #shell mesh
    if mesh.geometry().dim() == mesh.topology().dim() + 1:
        for v in vertices(mesh):

            dof = v2d[v.index()]

            avg_normal = np.zeros(MaxDim)

            counter = 0.0

            for c in cells(v):

                if FacetMarker != None:

                    Use = False

                    marker = FacetMarker[c.index()]

                    for k in Region:

                        if marker in k:

                            Use = True

                else:

                    Use = True

                if Use:

                    if mesh.cell_orientations()[c.index()] < 0.5:

                        Factor = 1.0

                    else:

                        Factor = -1.0

                    for i in range(MaxDim):

                        avg_normal[i] += Factor*c.cell_normal()[i]

                    counter = counter + 1.0

            if counter > 0.0:

                length = 0.0

                for i in range(MaxDim):

                    length += (avg_normal[i]/counter)*(avg_normal[i]/counter)

                length = sqrt(length)

                for i in range(MaxDim):

                    values[dof[i]] = avg_normal[i]/counter/length

    #full mesh

    else:

        for v in vertices(mesh):

            avg_normal = np.zeros(MaxDim)

            counter = 0.0

            dof = v2d[v.index()]

            for f in facets(v):

                if f.exterior() == True:

                    if FacetMarker != None:

                        Use = False

                        marker = FacetMarker[c.index()]

                        for k in Region:

                            if marker in k:

                                Use = True

                    else:

                        Use = True

                    if Use:

                        for i in range(MaxDim):

                            avg_normal[i] += f.normal()[i]

                        counter = counter + 1.0

            if counter > 0.0:

                length = 0.0

                for i in range(MaxDim):

                    length += (avg_normal[i]/counter)*(avg_normal[i]/counter)

                length = sqrt(length)

                for i in range(MaxDim):

                    values[dof[i]] = avg_normal[i]/counter/length

    N.vector().set_local(values)

    N.vector().apply("")

    return N

def regularize_solution(u, u_max=0, i_list = None):
    u_array = u.vector().get_local()
    
    if i_list is None:
        i_list = []

        for i,element in enumerate(u_array):
            if abs(element) > u_max: 
                u_array[i]=0
                i_list.append(i)
        u.vector()[:] = u_array
        return u, i_list
    else:
        for i in i_list:
            u_array[i]=0
            print(i)
        u.vector()[:] = u_array
        return u

def tang_grad(u, normals):
    return grad(u) - outer(dot(grad(u), normals), normals)

def get_objective(target_geometry, mesh, S0, d):
    eps = 1e-10
    if d == 2 and target_geometry == "circle":
        theta = Expression('atan2((x[1]),(x[0]))', degree = 1)
        return Expression(('S0*(cos(theta)*cos(theta)-0.5)', 'S0*sin(theta)*cos(theta)'), theta = theta, S0 = S0, degree = 1)

    elif d == 2 and target_geometry == "uniform_vertical":
        return Expression(('-S0*(0.5)', 'eps'), S0 = S0, eps=eps, degree = 1)
   

    elif d == 3 and target_geometry == "pseudoChiral":
        alpha = np.pi/4 # The angle of the director at the top of the domain
        h = 5 # The height of the domain
        return Expression(('S0*(cos(phi)*cos(phi)-1/3)', 'S0*sin(phi)*cos(phi)', 'eps', 'S0*sin(phi)*sin(phi) - 1/3', 'eps'), phi = Expression('alpha*x[2]/h', alpha=alpha, h=h, degree = 1) , S0 = S0, eps= eps, degree = 1)

    elif d == 3 and target_geometry == "uniform_horizontal":
        phi = Constant(np.pi/4)
        # Here we dont use Expression since it is less accurate than as_vector. Both are equivalent, since the target field is constant throughout space, and both implementations have been tested to give the same result up to numerical precision.
        return as_vector((S0*(cos(phi)*cos(phi)-1/3), S0*sin(phi)*cos(phi), eps, S0*sin(phi)*sin(phi) - 1/3, eps))
        # q_target = Expression(('S0*(cos(phi)*cos(phi)-1/3)', 'S0*sin(phi)*cos(phi)', 'eps', 'S0*(sin(phi)*sin(phi) - 1/3)', 'eps'), phi = phi , S0 = S0, eps= eps, degree = 1)



    else:
        raise ValueError("Target geometry not supported.")