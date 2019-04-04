from firedrake import *
from firedrake.petsc import PETSc
import numpy as np
import math


class Problem:
    def __init__(self,sedimentation_rate,fluid_source,initial_depth):
        ''' Initialize fluid and rock properties '''
        # General Terms
        self.y2s = 365.0*24.0*60.0*60.0 # Number of seconds in 1 year
        self.g = Constant((0.,0.,9.81)) # this represents g*z where z is the upward pointing unit vector        
        self.g_constant = 9.81 # Acceleration due to gravity (scalar)
        self.initial_depth = initial_depth # Initial depth used to compute overburden stress and initial hydrostatic fluid pressure for DirichletBC
        self.rho_ref = 2700. # Reference density for overburden stress given initial depth kg/m^3
        # Specify Fluid Properties
        self.mu = 0.0005 # Viscosity Pa-s (0.005 Poise or 0.5 cP), 1 Pa-s = 10 Poise
        self.rho_fluid = 1100. # Fluid Density
        # Specify Elastic Rock Properties
        self.phi_sand = 0.20
        self.rho_sand = 2650. # Density of sandstone = 2650 kg/m^3
        self.poisson = 0.25 # Poisson's ratio (denoted as nu in derivation)
        self.bulk_compressibility = 1.0e-8 # Bulk rock compressibility [Pa]
        self.lambda_elastic = self.bulk_compressibility*((1.0+self.poisson)/(3.0*(1.0-self.poisson)))
        self.k_sand = 2.46730816679E-13 # Sand permability = 250 mD -> m^2
        # Specify Viscous Rock Properties
        self.graindiameter = 0.001 # Grain diameter = 1 mm
        self.salt_thickness = 10.0
        self.rho_salt = 2060. # Density of rock salt = 2060 kg/m^3        
        self.phi_salt = 0.005
        self.viscosity_ref = 1.0e17 # Reference effective viscosity of the ductile layer (Pa-s) 
        self.m = 1 # Exponent for porosity in the ductile layer
        self.k_salt = ((self.phi_salt**2)*(self.graindiameter**2))/1600 #4.0e-20 # Salt permeability = 40 Nanodarcy-> m^2
        # Specify Source Terms
        self.rho_sed = 2700. # Density of sediment deposited, Sandstone = 2700 kg/m^3
        self.sed_rate = sedimentation_rate # Sedimentation rate (denoted as greek letter Psi in derivation)
        self.omega = self.rho_sed*self.sed_rate*self.g_constant # Contribution due to overlying sedimentation
        self.fs = fluid_source # Fluid source term (denoted as fs in derivation)


    def import_mesh(self,path_to_mesh_file,element_order):
        """ Imports the mesh from a .msh file. Element order must agree with
        that in the meshfile to establish appropriate function spaces"""
        PETSc.Sys.Print('Setting up mesh across %d processors' % COMM_WORLD.size)
        self.mesh = Mesh(path_to_mesh_file,comm=COMM_WORLD)
        PETSc.Sys.Print('Mesh imported from '+path_to_mesh_file)
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
        self.top_ss_thickness = 20.0 # Thickness of top Sandstone unit

    def initialize_functions(self):
        ''' Create trial, test functions as well as indicator functions for
        subdomains '''
        PETSc.Sys.Print('Initializing functions...')
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
        self.indicator_all =  Function(funcspace_constant,name="Indicator-all-Volume")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] =1.0;',dx,{'f':(self.indicator_all, WRITE)} )
        self.indicator_sandall = Function(funcspace_constant,name="Indicator-Sand-All(no interface)").interpolate(self.indicator_all-self.indicator_salt)
        self.indicator_sandtop_ni = Function(funcspace_constant,name="Indicator-Sand-Top(no interface)")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 1.0;',dx(self.top_sand_vol_id), {'f':(self.indicator_sandtop_ni, WRITE)} )
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 0.0;',dx(self.salt_vol_id), {'f': (self.indicator_sandtop_ni, WRITE)} )
        self.indicator_sandbottom_ni = Function(funcspace_constant,name="Indicator-Sand-Bottom(no interface)")
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 1.0;',dx(self.bottom_sand_vol_id), {'f':(self.indicator_sandbottom_ni, WRITE)} )
        par_loop( 'for ( int i=0; i < f.dofs; i++ ) f[i][0] = 0.0;',dx(self.salt_vol_id), {'f': (self.indicator_sandbottom_ni, WRITE)} )
        self.porosity = Function(FunctionSpace(self.mesh,"CG",3),name="Porosity").interpolate(self.indicator_sandall*self.phi_sand+self.indicator_salt*self.phi_salt)
        self.k = Function(FunctionSpace(self.mesh,"CG",3),name="Permeability").interpolate(self.indicator_sandall*self.k_sand+self.indicator_salt*self.k_salt)

        # Create function for overburden stress
        x,y,z = SpatialCoordinate(self.mesh) # Extract coordinates of mesh
        funcspace_discontinuous_linear = FunctionSpace(self.mesh,"DG",1) # create function space using linear piecewise polynomial
        funcspace_continuous_linear = FunctionSpace(self.mesh,"CG",1)
        reference_overburden = self.initial_depth*self.g_constant*self.rho_ref
        self.overburden = Function(funcspace_continuous_linear).interpolate((-self.rho_sand*self.g_constant*z+reference_overburden)*self.indicator_sandtop_ni + \
                                                                            (-self.rho_salt*self.g_constant*(z+self.top_ss_thickness) \
                                                                            +self.rho_sand*self.g_constant*self.top_ss_thickness+reference_overburden)*self.indicator_salt \
                                                                            +(self.salt_thickness*self.rho_salt*self.g_constant+self.rho_sand*self.g_constant*self.top_ss_thickness \
                                                                            + reference_overburden -self.rho_sand*self.g_constant*(z+self.top_ss_thickness+self.salt_thickness))*self.indicator_sandbottom_ni) # Overburden stress


    def boundary_conditons(self,dirichlet_value):
        ''' Defines boundary conditions for top and bottom boundary
        where dirichlet and neumann are tuples representing (value,index) '''
        PETSc.Sys.Print('Setting up boundary conditions')
        self.BCDirichlet = Constant(dirichlet_value) # Value of Dirichlet Boundary Condition
        self.BCDirichletIndex = self.top_id # Index of Dirichlet Boundary Condition
        self.BCNeumann = Constant(9.81*self.rho_fluid) # Value of Neumann Boundary Condition
        self.BCNeumannIndex = self.bottom_id # Index of Neumann Boundary Condition
        self.bcs = DirichletBC(self.funcspace,self.BCDirichlet,self.BCDirichletIndex) # Create Boundary Condition


    def update_parameters(self):
        ''' Updates porosity and permeability based on the evolution of the fluid pressure'''
        phi_salt = self.dt*(self.porosity*self.Pressure/self.viscosity_ref - self.porosity*self.overburden/(self.viscosity_ref*(1-self.porosity)) \
                            + self.porosity*self.porosity*self.Pressure/(self.viscosity_ref*(1-self.porosity))) + self.porosity
        phi_sand = self.phi_ref + self.bulk_compressibility*(1-self.phi_ref)*(self.Pressure - self.Pressure_ref)
        self.porosity = Function(FunctionSpace(self.mesh,"CG",3),name="Porosity").interpolate(self.indicator_sandall*phi_sand+self.indicator_salt*phi_salt)
        k_salt = ((phi_salt*phi_salt)*(self.graindiameter**2))/1600
        
    
    def solve_system(self):
        self.rho=self.rho_fluid
        PETSc.Sys.Print('Solving problem ...')
        self.k = self.k_sand
        a = dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.bottom_sand_vol_id) + \
        dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.salt_vol_id) + \
        dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.top_sand_vol_id)
        L = (self.k/self.mu)*self.BCNeumann*self.V*ds(self.BCNeumannIndex) - \
        dot((self.k/self.mu)*self.g*self.rho, grad(self.V))*dx \
        + (self.k/self.mu)*self.rho*dot(self.g,self.normal)*self.V*ds(self.BCNeumannIndex)
        self.Pressure = Function(self.funcspace)
        solve(a == L, self.Pressure, bcs=self.bcs, solver_parameters={"ksp_type": "gmres"})
        #F = a-L
        #solve(F == 0, self.Pressure, bcs = self.bcs) # For non-linear variation problems

    
    def solve_viscous(self):
        ''' Solves for the fluid pressure in a viscous domain
        Note: This solves the viscous equations over the entire domain'''
        PETSc.Sys.Print('Solving viscous problem...')
        self.porosity = 0.01
        self.k = self.k_salt
        self.viscosity_ref = 1.0e10
        a = (dot(grad(self.V),(self.k/self.mu)*grad(self.U))+((self.porosity/self.viscosity_ref)*self.U*self.V))*dx
        L =((((self.overburden*self.porosity)/(self.viscosity_ref*(1.0-self.porosity)))*self.V)-dot(grad(self.V),(self.k/self.mu)*self.rho_fluid*self.g))*dx \
                +(self.k/self.mu)*self.BCNeumann*self.V*ds(self.BCNeumannIndex) + (self.k/self.mu)*self.rho_fluid*dot(self.g,self.normal)*self.V*ds(self.BCNeumannIndex)
        self.PressureViscous = Function(self.funcspace,name="Fluid Pressure Viscous")
        solve(a == L, self.PressureViscous, bcs=self.bcs, solver_parameters={'ksp_type':\
                                                                    "gmres",'ksp_gmres_restart':200,"ksp_rtol":\
                                                                    1e-10,'ksp_monitor_singular_value':\
                                                                    None})
    def solve_elastic(self):
        ''' Solves for the fluid pressure (steady-state) in an elastic domain
        Note: This solves the elastic equations over the entire domain'''
        PETSc.Sys.Print('Solving elastic problem...')
        a = dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.bottom_sand_vol_id) + \
        dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.salt_vol_id) + \
        dot((self.k/self.mu)*grad(self.U), grad(self.V))*dx(self.top_sand_vol_id)
        L = (self.k/self.mu)*self.BCNeumann*self.V*ds(self.BCNeumannIndex) - \
        dot((self.k/self.mu)*self.g*self.rho_fluid, grad(self.V))*dx \
        + (self.k/self.mu)*self.rho_fluid*dot(self.g,self.normal)*self.V*ds(self.BCNeumannIndex)
        self.PressureElastic = Function(self.funcspace, name="Fluid Pressure Elastic")
        solve(a == L, self.PressureElastic, bcs=self.bcs, solver_parameters={'ksp_type':\
                                                                    "gmres",'ksp_gmres_restart':200,"ksp_rtol":\
                                                                    1e-10,'ksp_monitor_singular_value':\
                                                                    None})

    def solve_instantaneous(self):
        ''' Solves for the instantaneous fluid pressure in the coupled domain '''
        PETSc.Sys.Print('Solving instantaneous problem...')
        # Contributions from bottom sandstone, a = L.H.S. and l = R.H.S.
        #self.bc2 = DirichletBC(self.funcspace,self.BCDirichlet+50.0*self.rho_fluid*self.g_constant,self.BCNeumannIndex)
        a1 = dot(grad(self.V),(self.k/self.mu)*grad(self.U))*dx(self.bottom_sand_vol_id)
        l1 = self.V*((self.k/self.mu)*(self.BCNeumann+self.rho_fluid*dot(self.g,self.normal)))*ds(self.BCNeumannIndex)+ \
        (self.V*self.fs-(dot(grad(self.V),(self.k/self.mu)*self.rho_fluid*self.g)))*dx(self.bottom_sand_vol_id)
        # Contributions from middle salt layer
        a2 = (dot(grad(self.V),(self.k/self.mu)*grad(self.U))+((self.porosity/self.viscosity_ref)*self.U*self.V))*dx(self.salt_vol_id)
        l2=((((self.overburden*self.porosity)/(self.viscosity_ref*(1.0-self.porosity)))*self.V)-dot(grad(self.V),(self.k/self.mu)*self.rho_fluid*self.g))*dx(self.salt_vol_id)
        # Test to see if we can solve elastic
        #a2 = dot(grad(self.V),(self.k2/self.mu)*grad(self.U))*dx(self.salt_vol_id)
        #l2 = -dot(grad(self.V),(self.k2/self.mu)*self.rho_fluid*self.g)*dx(self.salt_vol_id)
        # Contributions from top sandstone
        a3 = dot(grad(self.V),(self.k/self.mu)*grad(self.U))*dx(self.top_sand_vol_id)
        l3 = -dot(grad(self.V),(self.k/self.mu)*self.rho_fluid*self.g)*dx(self.top_sand_vol_id)

        A = a1+a2+a3
        L = l1+l2+l3
        F = A-L
        self.PressureCoupled = Function(self.funcspace,name="Fluid Pressure Coupled")
        solve(A == L,self.PressureCoupled,bcs=[self.bcs],solver_parameters={'ksp_type':\
                                                                    "gmres",'ksp_gmres_restart':200,"ksp_rtol":\
                                                                    1e-12})
        #solve(F == 0, self.Pressure, bcs = self.bcs) # For non-linear variation problems


    def comp_flux_interface(self,Function):
        ''' Computes the flux across the salt bottom interface. Note this is
        reported as a volumetric flux and should be scaled by the porosity '''
        PETSc.Sys.Print('Computing flux across salt interface')
        self.flux_across_interface = assemble(dot((-self.k/self.mu)*(grad(Function)+self.rho_fluid*self.g),self.normal)('+')*dS(self.salt_bottom_id))
        return self.flux_across_interface

    def save_solution(self,filename,function):
        File(filename).write(function)
