from fenics import *



n = 30
eps = 0.01


#create the mesh*

mesh = Mesh('disc.xml')



#define vector spaces*
V = VectorFunctionSpace(mesh, 'Lagrange', 1, dim=2)
W = FunctionSpace(mesh,'Lagrange',1)


#define boundary conditions
f = Expression(('-1*x[1]','x[0]'),degree = 2)
def boundary(x,on_boundary):
    return on_boundary

bc = DirichletBC(V,f,boundary)

#define trial and test function*
u = TrialFunction(V)
v1, v2 = TestFunctions(V)
u1, u2  = split(u)

#define variational problem
a = eps*dot(grad(u1),grad(v1))*dx+dot(grad(u2),grad(v2))*dx+(u1*u1+u2*u2-1)*(u1*v1+u2*v2)/eps*dx


#solve variational problem
u_ = Function(V)
F = action(a,u_)

J = derivative(F,u_,u)

problem = NonlinearVariationalProblem(F,u_,bc,J)
solver = NonlinearVariationalSolver(problem)

solver.parameters['newton_solver']['absolute_tolerance']=1E-6
solver.parameters['newton_solver']['maximum_iterations']=10000

solver.solve()




#save the solution
u1, u2 = u_.split()
s = interpolate(Expression('pow(f1*f1+f2*f2,0.5)',f1=u1,f2=u2,degree=2),W)

file = File('u1.pvd') 
file << u1
file = File('u2.pvd')
file << u2
file = File('s.pvd')
file << s


#plot the solution

