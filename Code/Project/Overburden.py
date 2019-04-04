from TrapIntegrity import *

meshfile = '../Mesh/Box/MediumResolution/3DDomain.msh'
Example = Problem(sedimentation_rate=0.,fluid_source=0.,initial_depth = 10.)
Example.import_mesh(meshfile,element_order=2)
Example.initialize_functions()
print('Flux: ' + Example.comp_flux_interface(Example.Pressure))
Example.save_solution('Permeability.pvd',Example.k)
Example.save_solution('Porosity.pvd',Example.porosity)
