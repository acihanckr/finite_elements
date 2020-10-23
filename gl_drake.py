from firedrake import *

n = 30
eps = 0.01


#create the mesh*
mesh = UnitDiskMesh(n)


#define vector spaces*
V = VectorFunctionSpace(mesh, 'Lagrange', 1, dim=2)
W = FunctionSpace(mesh,'Lagrange',1)


#define boundary conditions
f = Expression(('-sin(x[1])','cos(x[0])'))
bc = DirichletBC(V,f,1)

#define trial and test function*
u = TrialFunction(V)
v1, v2 = TestFuncions(V)
u1, u2  = split(u)

#define variational problem
a = eps*dot(grad(u1),grad(v1))+dot(grad(u2),grad(v2))+(u1**2+u2**2-1)*u*v/eps*dx


#solve variational problem


#save the solution


#plot the solution
