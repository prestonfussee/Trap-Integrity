
# Import FEniCS
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import math

# Top Boundary: Helper function required by fenics
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2],0.0)

# Bottom Boundary: Helper function required by fenics
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2],-20.0)

# Left Boundary: Helper function required by fenics
class SideLeft(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],0.0)

# Right Boundary: Helper function required by fenics
class SideRight(SubDomain): 
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 20.0)

# Front Boundary: Helper function required by fenics
class SideFront(SubDomain): 
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1],0.0)

# Back Boundary: Helper function required by fenics
class SideBack(SubDomain): 
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1],20.0)
    
class FlowElastic3D:
    def __init__(self):
        self.d = 2 # order of the elements, here quadratic piecewise polynomials are used
        self.g = Constant((0.,0.,9.81)) # this represents g*z where z is the upward pointing unit vector
        self.k = 2.46730816679E-13 # 250 mD convert to m^2
        self.mu = 0.0005 # 0.5 cP (Viscosity value for a sodium formate brine w/ rho = 1100 and T = 100 C )
        self.rho = Constant(1100.) # Density of 1100 kg/m^3
        self.generate_mesh() # Generates the mesh
        self.set_boundary_conditions() # Enforces Boundary Conditions
        
    def generate_mesh(self):
        x0 = 0. # Define min(x)
        y0 = 0. # Define min(y)
        x1 = 20. # Define max(x)
        y1 = 20. # Define max(y)
        z1 = 0. # Define max(z)
        z0 = -20. # Define min(z)
        nx = 25 # Number of cells in x-direction
        ny = 25 # Number of cells in y-direction
        nz = 25 # Number of cells in z-direction
        self.mesh = BoxMesh(Point(x0, y0, z0), Point(x1, y1, z1), nx, ny, nz) # creates rectangular mesh
        self.Vh = FunctionSpace(self.mesh, "Lagrange", self.d) # Define function space on mesh (Lagrange with order = self.d)
        # Initialize mesh function for interior domains
        self.domains = MeshFunction('size_t', self.mesh, self.mesh.topology().dim()-1)
        self.domains.set_all(0) # Set all cell markets to 0 (Later changed in set_boundary_conditions)
    
        
    def set_boundary_conditions(self):
        # Instantiate classes that prescribe boundaries
        self.TopDirichlet  = Constant(0.) # Top most boundary condition (Dirichlet)
        self.SideNeumann = Constant(0.) # Boundary condition on sides (Neumann)
        self.BottomNeumann = Constant(9.81*1100.0) # Bottom boundary Condition = rho*g
        self.TopBoundary = Top() # Instantiate classes
        self.BottomBoundary = Bottom() # Instantiate classes
        self.LeftBoundary = SideLeft() # Instantiate classes
        self.RightBoundary = SideRight() # Instantiate classes
        self.FrontBoundary = SideFront() # Instantiate classes
        self.BackBoundary = SideBack() # Instantiate classes
        self.TopBoundary.mark(self.domains, 1) # Mark cells belonging to top boundary with integer = 1
        self.BottomBoundary.mark(self.domains, 2) # Mark cells belonging to bottom boundary with integer = 2
        self.LeftBoundary.mark(self.domains, 3) # Mark cells belonging to left side boundary with integer = 3
        self.RightBoundary.mark(self.domains, 4) # Mark cells belonging to right side boundary with integer = 4
        self.FrontBoundary.mark(self.domains, 5) # Mark cells belonging to front side boundary with integer = 5
        self.LeftBoundary.mark(self.domains, 6) # Mark cells belonging to back side boundary with integer = 6                       
        self.bcs = [DirichletBC(self.Vh, self.TopDirichlet, self.domains, 1)] # Create object to store Dirichlet Boundary Conditions
        self.ds = Measure("ds", subdomain_data=self.domains) # Create a measure to easily identify parts of the domain

    def solve_system(self):
        # Construct function spaces to search for solution
        n = FacetNormal(self.mesh) # Get normal component to mesh
        uh = TrialFunction(self.Vh) # Set up trial function based on function space (Vh)
        vh = TestFunction(self.Vh) # Set up test function based on funtion space(Vh)
        a = inner((self.k/self.mu)*grad(uh), grad(vh))*dx # Define left hand side 
        L = (self.k/self.mu)*self.BottomNeumann*vh*self.ds(2) - inner((self.k/self.mu)*self.g*self.rho, grad(vh))*dx \
        + (self.k/self.mu)*self.rho*dot(self.g,n)*vh*self.ds(2) # Define right hand side
        A, b = assemble_system(a, L, self.bcs) # Assembly
        self.u = Function(self.Vh) # Define place holder function
        solve(A, self.u.vector(), b) # Solve for weights for function u
        # Return the solution
        return self.u
      
        
    def plot_solution(self):
        # Define plot size
        plt.rcParams['figure.figsize'] = [22, 11]
        # Plot Numerical Solution
        s = plot(self.solve_system()/1e6)
        ax = plt.gca()
        ax.set_xlabel('X [M]')
        ax.set_ylabel('Y [M]')
        ax.set_zlabel('Depth [M]')
        #plt.zlabel("Depth [M]")
        plt.title("Approximation of Hydrostatic Fluid Pressure")
        plt.colorbar(s,label='Pressure [MPa]')
        # Show Plots
        plt.show(s)
        

TestProb = FlowElastic3D()
x = TestProb.solve_system()
print(len(x))