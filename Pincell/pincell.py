import openmc
import matplotlib.pyplot as plt
import numpy as np
#OPENMC_CROSS_SECTIONS = "/Users/Julian/nndc_hdf5/cross_sections.xml"


###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
batches = 100
inactive = 8
particles = 5000

###############################################################################
#                 Exporting to OpenMC materials.xml file
###############################################################################


# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(material_id=1, name='UO2 fuel at 2.4% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=2.4)
uo2.add_element('O', 2.)

helium = openmc.Material(material_id=2, name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)

zircaloy = openmc.Material(material_id=3, name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')

#Borated water:
borated_water = openmc.Material(material_id=4, name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 670*1e-6)
borated_water.add_element('H', 2.0)
borated_water.add_element('O', 1.0)
borated_water.add_s_alpha_beta('c_H_in_H2O')


#Heavy water:
heavy_water = openmc.Material(material_id=5, name='Heavy water')
heavy_water.set_density('g/cm3', 1.10452)
heavy_water.add_nuclide('H2', 2.)
heavy_water.add_element('O', 1.)
heavy_water.add_s_alpha_beta('c_D_in_D2O')


#Normal water:
normal_water = openmc.Material(material_id=6, name="normal_water")
normal_water.set_density('g/cm3',0.99701)
normal_water.add_nuclide('H1',2.0)
normal_water.add_element('O',1.0)
normal_water.add_s_alpha_beta('c_H_in_H2O')

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, helium, zircaloy, borated_water,
                heavy_water, normal_water])
materials_file.export_to_xml()

###############################################################################
#                 Exporting to OpenMC geometry.xml file
###############################################################################

# Instantiate ZCylinder surfaces
#fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, r=0.39218, name='Fuel OR')
#clad_ir = openmc.ZCylinder(surface_id=2, x0=0, y0=0, r=0.40005, name='Clad IR')
#clad_or = openmc.ZCylinder(surface_id=3, x0=0, y0=0, r=0.45720, name='Clad OR')
fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, r=0.22, name='Fuel OR')
clad_ir = openmc.ZCylinder(surface_id=2, x0=0, y0=0, r=0.23, name='Clad IR')
clad_or = openmc.ZCylinder(surface_id=3, x0=0, y0=0, r=0.26, name='Clad OR')



left = openmc.XPlane(surface_id=4, x0=-0.62992, name='left')
right = openmc.XPlane(surface_id=5, x0=0.62992, name='right')
bottom = openmc.YPlane(surface_id=6, y0=-0.62992, name='bottom')
top = openmc.YPlane(surface_id=7, y0=0.62992, name='top')

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
top.boundary_type = 'reflective'
bottom.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(cell_id=1, name='cell 1')
gap = openmc.Cell(cell_id=2, name='cell 2')
clad = openmc.Cell(cell_id=3, name='cell 3')
water = openmc.Cell(cell_id=4, name='cell 4')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
gap.region = +fuel_or & -clad_ir
clad.region = +clad_ir & -clad_or
water.region = +clad_or & +left & -right & +bottom & -top

# Register Materials with Cells
fuel.fill = uo2
gap.fill = helium
clad.fill = zircaloy
water.fill = borated_water
#water.fill = heavy_water
#water.fill = normal_water

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
root.add_cells([fuel, gap, clad, water])

# Instantiate a Geometry, register the root Universe, and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()

###############################################################################
#                   Exporting to OpenMC settings.xml file
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {'tallies': True}

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 0.62992, 0.62992, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left = [-0.39218, -0.39218, -1.e50]
entropy_mesh.upper_right = [0.39218, 0.39218, 1.e50]
entropy_mesh.dimension = [10, 10, 1]
settings_file.entropy_mesh = entropy_mesh
settings_file.export_to_xml()

###############################################################################
#                   Exporting to OpenMC tallies.xml file
###############################################################################

# Instantiate a tally mesh
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100, 1]
mesh.lower_left = [-0.62992, -0.62992, -1.e50]
mesh.upper_right = [0.62992, 0.62992, 1.e50]

# Instantiate some tally Filters
energy_filter = openmc.EnergyFilter([0., 0.025,0.125,4., 10.,20.,100.,1000.,
                    10000.,50000.,20000000.])
mesh_filter = openmc.MeshFilter(mesh)
meshsurface_filter = openmc.MeshSurfaceFilter(mesh)

particle_filter = openmc.ParticleFilter(['neutron','photon'])

# Instantiate the Tally
tally = openmc.Tally(tally_id=1, name='tally 1')
tally.filters = [energy_filter, mesh_filter, particle_filter]
tally.scores = ['flux', 'fission', 'nu-fission','prompt-nu-fission']
'''
#tally for neutron flux distribution plot
energy_bins2 = np.logspace(-3,7,50)
energyfilter2 = openmc.EnergyFilter(energy_bins2)
tally_energy = openmc.Tally(tally_id=99, name='energy flux tally')
tally_energy.filters = [energyfilter2]#,particle_filter]
tally_energy.scores = ['flux']
'''

# Instantiate a Tallies collection and export to XML
tallies_file = openmc.Tallies([tally])#, tally_energy])

###############################################################################
#                            Four factor formula terms                        #
###############################################################################

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.025]),
                                openmc.CellFilter([fuel])]
tallies_file.append(fuel_therm_abs_rate)


# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.025])]
tallies_file.append(therm_abs_rate)


# Fast Fission Factor tallies
therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
therm_fiss_rate.scores = ['nu-fission']
therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.025])]
tallies_file.append(therm_fiss_rate)

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
abs_rate = openmc.Tally(name='abs. rate')
fiss_rate.scores = ['nu-fission']
abs_rate.scores = ['absorption']
tallies_file += (fiss_rate, abs_rate)


# Instantiate thermal, fast, and total leakage tallies
leak = openmc.Tally(name='leakage')
leak.filters = [meshsurface_filter]
leak.scores = ['current']
tallies_file.append(leak)

thermal_leak = openmc.Tally(name='thermal leakage')
thermal_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0., 0.025])]
thermal_leak.scores = ['current']
tallies_file.append(thermal_leak)

fast_leak = openmc.Tally(name='fast leakage')
fast_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0.025, 20.0e6])]
fast_leak.scores = ['current']
tallies_file.append(fast_leak)


# Instantiate a Tallies collection and export to XML
tallies_file.export_to_xml()


p = openmc.Plot()
p.filename = 'materials-xy'
p.origin = [0, 0, 0]
p.basis = 'xy'
p.width = [1.0, 1.0]
p.pixels = [3000, 3000]
p.color_by = 'material'
p.colors = {heavy_water:'cyan',zircaloy: 'burlywood', helium: 'limegreen',
            uo2:'yellow'}

plots = openmc.Plots([p])
plots.export_to_xml()

openmc.run()
#openmc.plot_geometry()


###############################################################################
#                            Tally data processing                            #
###############################################################################

sp = openmc.StatePoint('statepoint.100.h5')
sp_tally_flux = sp.get_tally(scores=['flux'])
sp_tally_fission = sp.get_tally(scores=['fission'])

flux = sp_tally_flux.get_slice(scores=['flux'])
thermal_flux = flux.get_slice(filters=[openmc.EnergyFilter],
                filter_bins=[((0.,0.025),)])
thermal_neutron_flux = thermal_flux.get_slice(filters=[openmc.ParticleFilter],
                filter_bins = [(('neutron'),)], squeeze = True)


neutron_flux = flux.get_slice(filters=[openmc.ParticleFilter,
                openmc.EnergyFilter], filter_bins = [(('neutron'),),
                ((50000.,20000000.),)], squeeze = True)


fission = sp_tally_fission.get_slice(scores=['fission'])
therm_fission = fission.get_slice(filters=[openmc.EnergyFilter],
                filter_bins=[((0.,0.025),)])
therm_neutron_fission = therm_fission.get_slice(filters=[openmc.ParticleFilter],
                filter_bins = [(('neutron'),)], squeeze = True)


thermal_neutron_flux.std_dev.shape=(100,100)
thermal_neutron_flux.mean.shape=(100,100)

neutron_flux.std_dev.shape=(100,100)
neutron_flux.mean.shape=(100,100)

therm_neutron_fission.std_dev.shape=(100,100)
therm_neutron_fission.mean.shape=(100,100)


fig = plt.subplot(121)
fig.imshow(thermal_neutron_flux.mean)
fig2 = plt.subplot(122)
fig2.imshow(therm_neutron_fission.mean)
plt.show()

plt.imshow(neutron_flux.mean)
plt.colorbar()
plt.xlabel('x-direction')
plt.ylabel('y-direction')
#plt.title('Prompt neutron flux distribution')
#plt.savefig('promptfluxdist.png')
plt.show()


relative_error = np.zeros_like(thermal_flux.std_dev)
nonzero = thermal_flux.mean > 0
relative_error[nonzero]=thermal_flux.std_dev[nonzero]/(thermal_flux.mean[nonzero])
ret = plt.hist(relative_error[nonzero],bins=50)
plt.show()


'''tally_energy_data = sp.get_tally(name='energy flux tally')
flux_energy = tally_energy_data.get_slice(scores=['flux'])
bins = np.logspace(-3,7,49)
plt.loglog(bins,flux_energy.mean[:,0,0])
plt.xlabel('Neutron energy [eV]')
plt.ylabel('Flux [n/cm per source particle]')
plt.title(r'Energy distribution of the neutron flux using $H_{2}O$')
#plt.savefig('energydistfluxh2o.png',dpi=500)
plt.show()
'''

###############################################################################
#                            Source sites                                     #
###############################################################################

energy_bins = np.logspace(3,7)
probability, bin_edges = np.histogram(sp.source['E'], energy_bins,
    density=True)

print(sum(probability*np.diff(energy_bins)))

plt.semilogx(energy_bins[:-1], probability*np.diff(energy_bins),
    drawstyle='steps')
plt.xlabel('Energy [eV]')
plt.ylabel('Probability/eV')
plt.show()


plt.quiver(sp.source['r']['x'], sp.source['r']['y'],sp.source['u']['x'],
    sp.source['u']['y'], np.log(sp.source['E']), cmap='jet',scale=20.0)
plt.colorbar()
plt.xlim((-0.5,0.5))
plt.ylim((-0.5,0.5))
plt.show()

###############################################################################
#                           Printing Four factor formula terms                #
###############################################################################
# Get the fission and absorption rate tallies
fiss_rate = sp.get_tally(name='fiss. rate')
abs_rate = sp.get_tally(name='abs. rate')

# Get the leakage tally
leak = sp.get_tally(name='leakage')
leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

# Compute k-infinity using tally arithmetic
keff = fiss_rate / (abs_rate + leak)
print('k-infinity')
print(keff.get_pandas_dataframe())

# Compute resonance escape probability using tally arithmetic
therm_abs_rate = sp.get_tally(name='therm. abs. rate')
thermal_leak = sp.get_tally(name='thermal leakage')
thermal_leak = thermal_leak.summation(filter_type=openmc.MeshSurfaceFilter,
                remove_filter=True)
res_esc = (therm_abs_rate + thermal_leak) / (abs_rate + thermal_leak)
print('resonance escape')
print(res_esc.get_pandas_dataframe())

# Compute fast fission factor factor using tally arithmetic
therm_fiss_rate = sp.get_tally(name='therm. fiss. rate')
fast_fiss = fiss_rate / therm_fiss_rate
print('fast fission factor')
print(fast_fiss.get_pandas_dataframe())

therm_abs_rate = sp.get_tally(name='therm. abs. rate')
fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
therm_util = fuel_therm_abs_rate / therm_abs_rate
print('thermal Utilization')
print(therm_util.get_pandas_dataframe())

# Compute neutrons produced per absorption (eta) using tally arithmetic
eta = therm_fiss_rate / fuel_therm_abs_rate
print('eta')
print(eta.get_pandas_dataframe())

#Non-leakage terms
p_fnl = (abs_rate + thermal_leak) / (abs_rate + leak)
print('non-leakage fast')
print(p_fnl.get_pandas_dataframe())


p_tnl = therm_abs_rate / (therm_abs_rate + thermal_leak)
print('non-leakage thermal')
print(p_tnl.get_pandas_dataframe())

#K_eff
keff = res_esc * fast_fiss * therm_util * eta * p_fnl * p_tnl
print('full sixfactor k_eff')
print(keff.get_pandas_dataframe())
