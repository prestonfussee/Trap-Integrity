from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
import math


class Problem:
    def __init__(self):
        ''' Initialize fluid and rock properties '''
        self.g = Constant((0.,0.,9.81)) # this represents g*z where z is the upward pointing unit vector        
        self.k = 2.46730816679E-13 # Permeability m^2
        self.mu = 0.0005 # Viscosity Pa-s (0.005 Poise or 0.5 cP), 1 Pa-s = 10 Poise
        self.rho = 1100. # Fluid Density
    
    #@profile
    def import_mesh(self,path_to_mesh_file,element_order):
        """ Imports the mesh from a .msh file. Element order must agree with
        that in the meshfile to establish appropriate function spaces"""
        PETSc.Sys.Print('Setting up mesh across %d processors' % COMM_WORLD.size)
        self.mesh = Mesh(path_to_mesh_file,comm=COMM_WORLD)
        self.ele_order = element_order
        self.normal = FacetNormal(self.mesh) # normal component to mesh
        # Physical Groups
        # These ID numbers must match the ones defined in the mesh file
        self.bottom_id = 1 # Bottom Boundary
        self.salt_bottom_id = 2 # Salt Bottom Interface
        self.salt_top_id = 3  # Salt Top Interface
        self.top_id = 4 #Top Boundary
        self.side_id = 5 # Side Boundaries (all)
        self.bottom_sand_vol_id = 6 # Bottom Sandstone Volume
        self.salt_vol_id = 7  # Salt Volume
        self.top_sand_vol_id = 8 # Top Sandstone Volume
        self.funcspace = FunctionSpace(self.mesh,"CG",element_order) # Define function space

    #@profile
    def initialize_functions(self):
        ''' Create trial, test functions as well as indicator functions for
        subdomains '''
        self.U = TrialFunction(self.funcspace)
        self.V = TestFunction(self.funcspace)
        ### Create indicator functions that are 0 in one subdomain and 1 in the other
        funcspace_constant = FunctionSpace(self.mesh,"DG",0)
        self.indicator_sandbottom = Function(funcspace_constant,name="Indicator-Bottom-Sandstone-Volume")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 1.0;', dx(self.bottom_sand_vol_id), {'f': (self.indicator_sandbottom, WRITE)} )
        self.indicator_sandtop = Function(funcspace_constant,name="Indicator-Top-Sandstone-Volume")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 1.0;', dx(self.top_sand_vol_id), {'f': (self.indicator_sandtop, WRITE)} )
        self.indicator_salt = Function(funcspace_constant,name="Indicator-Salt-Volume")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 1.0;',dx(self.salt_vol_id), {'f': (self.indicator_salt, WRITE)} )
        ### Save Indicator Functions for Viewing
        #File('IndicatorSalt.pvd').write(self.indicator_salt)
        #File('IndicatorSandBottom.pvd').write(self.indicator_sandbottom)
        #File('IndicatorSandTop.pvd').write(self.indicator_sandtop)

    #@profile
    def boundary_conditons(self,dirichlet_value):
        ''' Defines boundary conditions for top and bottom boundary
        where dirichlet and neumann are tuples representing (value,index) '''
        PETSc.Sys.Print('Setting up boundary conditions')
        self.BCDirichlet = Constant(dirichlet_value) # Value of Dirichlet Boundary Condition
        self.BCDirichletIndex = self.top_id # Index of Dirichlet Boundary Condition
        self.BCNeumann = Constant(9.81*self.rho) # Value of Neumann Boundary Condition
        self.BCNeumannIndex = self.bottom_id # Index of Neumann Boundary Condition
        self.bcs = DirichletBC(self.funcspace,self.BCDirichlet,self.BCDirichletIndex) # Create Boundary Condition
    
    #@profile
    def solve_system(self):
        PETSc.Sys.Print('Solving problem ...')
        a1 = dot((self.k/self.mu)*self.indicator_sandbottom*grad(self.U),grad(self.V))*dx(self.bottom_sand_vol_id)
        a2 = dot((self.k/self.mu)*self.indicator_salt*grad(self.U),grad(self.V))*dx(self.salt_vol_id)
        a3 = dot((self.k/self.mu)*self.indicator_sandtop*grad(self.U),grad(self.V))*dx(self.top_sand_vol_id)
        l1 = (self.k/self.mu)*self.BCNeumann*self.V*ds(self.BCNeumannIndex) + (self.k/self.mu)*self.rho*dot(self.g,self.normal)*self.V*ds(self.BCNeumannIndex)
        l2 = -dot((self.k/self.mu)*self.g*self.rho, grad(self.V))*dx
        A = a1+a2+a3
        L = l1+l2
        self.Pressure = Function(self.funcspace)
        solve(A == L, self.Pressure, bcs=self.bcs, solver_parameters={"ksp_type": "gmres"})
        #F = a-L
        #solve(F == 0, self.Pressure, bcs = self.bcs) # For non-linear variation problems

    #@profile
    def comp_flux_interface(self):
        ''' Computes the flux across the salt bottom interface. Note this is
        reported as a volumetric flux and should be scaled by the porosity '''
        PETSc.Sys.Print('Computing flux across salt interface')
        self.flux_across_interface = -assemble(dot(grad(self.Pressure),self.normal)('+')*dS(self.salt_bottom_id))
        return self.flux_across_interface

    def save_solution(self,filename):
        File(filename).write(self.Pressure)

meshfile = '../Mesh/Box/MediumResolution/3DDomain.msh'
outputfile = 'HydrostaticPressurev2.pvd'
P = Problem()
P.import_mesh(meshfile,element_order=2)
P.initialize_functions()
P.boundary_conditons(0.0)
P.solve_system()
x = P.comp_flux_interface()
#PETSc.Sys.Print('Flux across salt interface = %g' % (P.comp_flux_interface()))
P.save_solution(outputfile)
