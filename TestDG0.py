from fenics import *
from dolfin import *
from dolfin_adjoint import *
from ufl_legacy import nabla_div, nabla_grad, VectorElement, FiniteElement, MixedElement, split, atan_2, replace, Jacobian

mesh = UnitSquareMesh(4, 4)


W = FunctionSpace(mesh, "DG", 1)
q = Expression('x[0]*x[1]', degree=1)
q = interpolate(q, W)
File("q.pvd") << q

gradq = grad(q)
proj_gradq = project(gradq, VectorFunctionSpace(mesh, "DG", 1))
for i, dof in enumerate(proj_gradq.vector().get_local()):
    print(f"DOF {i}: {dof}")
# proj_gradq = project(gradq, VectorFunctionSpace(mesh, "DG", 1))
# File("gradq.pvd") << proj_gradq