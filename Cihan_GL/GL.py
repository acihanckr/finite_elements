### Vector order parameter  ###

from dolfin import*

parameters["allow_extrapolation"]=True

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# Generate mesh
mesh2d = Mesh("GL.xml")

# Define function space
V = VectorFunctionSpace(mesh2d, "CG",1,dim=2)
V2 = FunctionSpace(mesh2d, "CG",1)

# Boundary conditions
gint = Expression(('x[0]*pow(x[0]*x[0]+x[1]*x[1],-0.5)','x[1]*pow(x[0]*x[0]+x[1]*x[1],-0.5)'),degree=2)
gext = Expression(('1.0','0.0'),degree=2)
def GammaI(x, on_boundary):
	r = pow(x[0] * x[0] + x[1] * x[1],0.5)
	return near(r,1.0,0.01)
def GammaE(x, on_boundary):
	r = pow(x[0] * x[0] + x[1] * x[1],0.5)
	return near(r,5.0,0.1)
bcint = DirichletBC(V,gint,GammaI)
bcext = DirichletBC(V,gext,GammaE)
bc = [bcint, bcext]

# Define variational problem
u=TrialFunction(V)
u1, u2 = split(u)
(v1, v2) = TestFunctions(V)
f = Constant("0.0")
a = inner(grad(u1), grad(v1))*dx + inner(grad(u2), grad(v2))*dx + 100.0*(u1*u1+u2*u2-1.0)*u1*v1*dx + 100.0*(u1*u1+u2*u2-1.0)*u2*v2*dx
L = f*v1*dx 

# Compute solution
#u_ = interpolate(Expression(('x[0]*x[0]*pow(x[0]*x[0]+x[1]*x[1],-1.0)-1.0/3.0','-1.0/3.0','0.0','x[0]*x[1]*pow(x[0]*x[0]+x[1]*x[1],-1.0)','0.0')),V) # initial guess
u_ = interpolate(Expression(('0.0','0.0'),degree=2),V) # initial guess
a = action(a, u_)
J = derivative(a, u_, u)

problem = NonlinearVariationalProblem(a, u_, bc, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-6
prm['newton_solver']['maximum_iterations'] = 10000

solver.solve()

u1, u2 = u_.split()
s = interpolate(Expression('pow(f1*f1+f2*f2,0.5)',f1=u1,f2=u2,degree=2),V2)

# Save solution
file = File("u1.pvd")
file << u1
file = File("u2.pvd")
file << u2
file = File("s.pvd")
file << s

