from TrapIntegrity import *
import math

meshfile = '../Mesh/Box/HighResolution/3DDomain.msh'
Salt_height = 10.0
Example = Problem(sedimentation_rate=0.,fluid_source=0.,initial_depth = 1000.)
Example.import_mesh(meshfile,element_order=2)
Example.boundary_conditons(1000.0*9.81*1100)
Total_depth = 1050.0
Hydrostatic = Total_depth*9.81*Example.rho_fluid
viscosity_ref = [1.0e10,1.0e11,1.0e12,1.0e13,1.0e14,1.0e15,1.0e16,1.0e17,1.0e18,1.0e19]
Porosity_values = np.linspace(0.001,0.20,100)
K = ((0.001*0.001)/1600)*Porosity_values*Porosity_values
Compaction_length = lambda k, viscosity: math.sqrt(k*viscosity/Example.mu)
Dimensionless_height = np.empty(len(Porosity_values)*len(viscosity_ref))
DimensionlessOverpressure = np.empty(len(Porosity_values)*len(viscosity_ref))
Fraction_of_verticalstress = np.empty(len(Porosity_values)*len(viscosity_ref))
counter = 0
for i in range(len(viscosity_ref)):
    for j in range(len(K)):
        print('Iteration: ' +str(counter)+'\n')
        Example.viscosity_ref = viscosity_ref[i]
        Example.phi_salt = Porosity_values[j]
        Example.k_salt = K[j]
        Example.initialize_functions()
        Example.solve_instantaneous()
        Dimensionless_height[counter] = Salt_height/Compaction_length(K[j],viscosity_ref[i])
        DimensionlessOverpressure[counter] = (Example.PressureCoupled((50,50,-50))-Hydrostatic)/Hydrostatic
        Fraction_of_verticalstress[counter] = Example.PressureCoupled((50,50,-50))/Example.overburden((50,50,-50))
        counter +=1

np.savetxt("DimensionlessHeight.csv", Dimensionless_height, delimiter=",")
np.savetxt("DimensionlessOverpressure.csv", DimensionlessOverpressure, delimiter=",")
np.savetxt("FractionofVerticalStress.csv",Fraction_of_verticalstress,delimiter=",")
#Example.save_solution('./Functions/Marc/FluidPressureNeumann.pvd',Example.PressureCoupled)
#Example.save_solution('./Functions/Marc/Lithostatic.pvd',Example.overburden#Example.save_solution('./Functions/FluidPressure.pvd',Example.Pressure)
#Example.save_solution('./Functions/OverburdenStress.pvd',Example.overburden)
#Example.save_solution('./Functions/Porosity.pvd',Example.porosity)
#Example.save_solution('./Functions/Permeability.pvd',Example.k)



