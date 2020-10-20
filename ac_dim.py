from __future__ import print_function
from fenics import *

eps = 0.01

mesh = UnitSquareMesh(8,8)
V = VectorFunctionSpace(mesh, 'P', 1,dim = 2)


def boundary(x, on_boundary):
    return on_boundary
#bc = DirichletBC(V, u_D, boundary)

u = TrialFunction(V)
u1, u2 = split(u)

v1,v2 = TestFunctions(V)

f = Constant(0.0)
a = grad(u1)[0]*grad(v1)[0]+grad(u2)[0]*grad(v2)[0]+grad(u1)[1]*grad(v1)[1]/(eps**2)+grad(u2)[1]*grad(v2)[1]/(eps**2)+grad(u1)[0]*grad(v1)[0]/eps+grad(u2)[1]*grad(v2)[1]/(eps**3)+4*(u1**2+u2**2-1)*(u1*v1+u2*v2)/(eps**2)*dx
L = f*v*dx

u = VectorFunction(V)
solve(a == L,u)


plot(u)
plot(mesh)

vtkfile = File('poisson/poisson/solution.pvd')
vtkfile<<u


error_L2 = errornorm(u_D, u, 'L2')


vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)

import numpy as np
error_max = np.max(np.abs(vertex_values_u_D-vertex_values_u))

print('error_L2 = ',error_L2)
print('error_max = ', error_max)

