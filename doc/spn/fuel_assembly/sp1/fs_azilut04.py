###############################################################################
## fs.py
## 9te [angband.ornl.gov]
## Wed Jan 12 10:37:50 2011
###############################################################################
## Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
##---------------------------------------------------------------------------##
## generated by /data/denovo/production/head/setup/bin/pygen built on 20110112
###############################################################################

import os, sys, math, string

# pyspn equation type
from spn_fv import *

print_it = False

##---------------------------------------------------------------------------##
## MAIN
##---------------------------------------------------------------------------##

initialize(sys.argv)

if node() == 0:
    print "Denovo - pyspn Python Front-End"
    print "-------------------------------"
    print "Release      : %16s" % (release())
    print "Release Date : %16s" % (release_date())
    print "Build Date   : %16s" % (build_date())
    print

timer = Timer()
timer.start()

##---------------------------------------------------------------------------##
## XS DATA

####### UO2 Fuel-Clad Macroscopic Cross Sections ##########
## Transport-corrected Total Cross Sections
T_UO2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
T_UO2[0] = 1.77949e-1
T_UO2[1] = 3.29805e-1
T_UO2[2] = 4.80388e-1
T_UO2[3] = 5.54367e-1
T_UO2[4] = 3.11801e-1
T_UO2[5] = 3.95168e-1
T_UO2[6] = 5.64406e-1

## Fission Cross Section
F_UO2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
F_UO2[0] = 7.21206e-3
F_UO2[1] = 8.19301e-4
F_UO2[2] = 6.45320e-3
F_UO2[3] = 1.85648e-2
F_UO2[4] = 1.78084e-2
F_UO2[5] = 8.30348e-2
F_UO2[6] = 2.16004e-1

## Nu
N_UO2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
N_UO2[0] = 2.78145
N_UO2[1] = 2.47443
N_UO2[2] = 2.43383
N_UO2[3] = 2.43380
N_UO2[4] = 2.43380
N_UO2[5] = 2.43380
N_UO2[6] = 2.43380

## Chi
C_UO2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
C_UO2[0] = 5.87910e-1
C_UO2[1] = 4.11760e-1
C_UO2[2] = 3.39060e-4
C_UO2[3] = 1.17610e-7
C_UO2[4] = 0.00000000
C_UO2[5] = 0.00000000
C_UO2[6] = 0.00000000

## Scattering Matrix for UO2 Fuel-Clad (Macroscopic)
S_UO2 = [ [[]], [[]], [[]], [[]], [[]], [[]], [[]]]
S_UO2[0] = [[1.27537e-1]]
S_UO2[1] = [[4.23780e-2], [3.24456e-1]]
S_UO2[2] = [[9.43740e-6], [1.63140e-3], [4.50940e-1]]
S_UO2[3] = [[5.51630e-9], [3.14270e-9], [2.67920e-3], [4.52565e-1], [1.25250e-4]]
S_UO2[4] = [[0.00000000], [0.00000000], [0.00000000], [5.56640e-3], [2.71401e-1], [1.29680e-3]]
S_UO2[5] = [[0.00000000], [0.00000000], [0.00000000], [0.00000000], [1.02550e-2], [2.65802e-1], [8.54580e-3]]
S_UO2[6] = [[0.00000000], [0.00000000], [0.00000000], [0.00000000], [1.00210e-8], [1.68090e-2], [2.73080e-1]]

## Upscattering Matrix
U_UO2 = [ [], [], [], [], [], [], [] ]
U_UO2[0] = []
U_UO2[1] = []
U_UO2[2] = []
U_UO2[3] = [4]
U_UO2[4] = [5]
U_UO2[5] = [6]
U_UO2[6] = []

######## 4.3% MOX Fuel-Clad Macroscopic Cross-Sections ############
## Transport-corrected Total Cross Sections
T_MOX43 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
T_MOX43[0] = 1.78731e-1
T_MOX43[1] = 3.30849e-1
T_MOX43[2] = 4.83772e-1
T_MOX43[3] = 5.66922e-1
T_MOX43[4] = 4.26227e-1
T_MOX43[5] = 6.78997e-1
T_MOX43[6] = 6.82852e-1

## Fission Cross-Sections
F_MOX43 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
F_MOX43[0] = 7.62704e-3
F_MOX43[1] = 8.76898e-4
F_MOX43[2] = 5.69835e-3
F_MOX43[3] = 2.28872e-2
F_MOX43[4] = 1.07635e-2
F_MOX43[5] = 2.32757e-1
F_MOX43[6] = 2.48968e-1

## Nu Cross-Sections
N_MOX43 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
N_MOX43[0] = 2.85209
N_MOX43[1] = 2.89099
N_MOX43[2] = 2.85486
N_MOX43[3] = 2.86073
N_MOX43[4] = 2.85447
N_MOX43[5] = 2.86415
N_MOX43[6] = 2.86780

## Chi Cross-Sections
C_MOX43 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
C_MOX43[0] = 5.87910e-1
C_MOX43[1] = 4.11760e-1
C_MOX43[2] = 3.39060e-4
C_MOX43[3] = 1.17610e-7
C_MOX43[4] = 0.00000000
C_MOX43[5] = 0.00000000
C_MOX43[6] = 0.00000000

## Scattering Matrix for 4.3% MOX Fuel-Clad (Macroscopic)
S_MOX43 = [ [[]], [[]], [[]], [[]], [[]], [[]], [[]] ]
S_MOX43[0] = [[1.28876e-1]]
S_MOX43[1] = [[4.14130e-2], [3.25452e-1]]
S_MOX43[2] = [[8.22900e-6], [1.63950e-3], [4.53188e-1]]
S_MOX43[3] = [[5.04050e-9], [1.59820e-9], [2.61420e-3], [4.57173e-1], [1.60460e-4]]
S_MOX43[4] = [[0.00000000], [0.00000000], [0.00000000], [5.53940e-3], [2.76814e-1], [2.00510e-3]]
S_MOX43[5] = [[0.00000000], [0.00000000], [0.00000000], [0.00000000], [9.31270e-3], [2.52962e-1], [8.49480e-3]]
S_MOX43[6] = [[0.00000000], [0.00000000], [0.00000000], [0.00000000], [9.16560e-9], [1.48500e-2], [2.65007e-1]]

## Upscattering Matrix
U_MOX43 = [ [], [], [], [], [], [], [] ]
U_MOX43[0] = []
U_MOX43[1] = []
U_MOX43[2] = []
U_MOX43[3] = [4]
U_MOX43[4] = [5]
U_MOX43[5] = [6]
U_MOX43[6] = []

############### Moderator 1 Macroscopic Cross-Sections ################
## Transport-corrected Total Cross Section
T_MOD1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
T_MOD1[0] = 1.59206e-1
T_MOD1[1] = 4.12970e-1
T_MOD1[2] = 5.90310e-1
T_MOD1[3] = 5.84350e-1
T_MOD1[4] = 7.18000e-1
T_MOD1[5] = 1.25445
T_MOD1[6] = 2.65038

## Scattering Matrix for Moderator (Macroscopic)
S_MOD1 = [ [[]], [[]], [[]], [[]], [[]], [[]], [[]] ]
S_MOD1[0] = [[4.44777e-2]]
S_MOD1[1] = [[1.13400e-1], [2.82334e-1]]
S_MOD1[2] = [[7.23470e-4], [1.29940e-1], [3.45256e-1]]
S_MOD1[3] = [[3.74990e-6], [6.23400e-4], [2.24570e-1], [9.10284e-2], [7.14370e-5]]
S_MOD1[4] = [[5.31840e-8], [4.80020e-5], [1.69990e-2], [4.15510e-1], [1.39138e-1], [2.21570e-3]]
S_MOD1[5] = [[0.00000000], [7.44860e-6], [2.64430e-3], [6.37320e-2], [5.11820e-1], [6.99913e-1], [1.32440e-1]]
S_MOD1[6] = [[0.00000000], [1.04550e-6], [5.03440e-4], [1.21390e-2], [6.12290e-2], [5.37320e-1], [2.48070   ]]

## Upscattering Matrix
U_MOD1 = [ [], [], [], [], [], [], [] ]
U_MOD1[0] = []
U_MOD1[1] = []
U_MOD1[2] = []
U_MOD1[3] = [4]
U_MOD1[4] = [5]
U_MOD1[5] = [6]
U_MOD1[6] = []

################### Create nuf vectors
NUF_UO2 = []
NUF_MOX43 = []
for i in range(0, 7):
   NUF_UO2.append( N_UO2[i] * F_UO2[i] )
   NUF_MOX43.append( N_MOX43[i] * F_MOX43[i] )

##---------------------------------------------------------------------------##
## BUILD MESH

def build_mesh(N):

    # vacuum    = 0
    # UO2       = 1
    # MOX       = 2
    # moderator = 3

    # UO2 pins
    uo2_pin  = Pincell()
    uo2_ids  = [1]
    uo2_r    = [0.4759]
    uo2_pin.set_shells(uo2_ids, uo2_r, 3)

    # MOX pins
    mox_pin  = Pincell()
    mox_ids  = [2]
    mox_r    = [0.4759]
    mox_pin.set_shells(mox_ids, mox_r, 3)

    # Make a 2x2 uo2 lattice and a 2x2 mox lattice
    uo2_lat = Lattice(2)
    mox_lat = Lattice(2)
    
    # lattices are uniform
    layout = [0, 0, 0, 0]
    uo2_lat.set_pins(layout)
    mox_lat.set_pins(layout)

    # assign the pins in the lattices
    uo2_lat.assign_pin(uo2_pin, 0)
    mox_lat.assign_pin(mox_pin, 0)

    # build the lattice
    uo2_lat.build_lattice(N)
    mox_lat.build_lattice(N)

    # print out mixing tables
    if print_it:
        print "UO2 Lattice"
        for m in xrange(uo2_lat.num_mixtures()):
            vf = uo2_lat.f(m)
            print "%4i" % (m),
            for f in vf:
                print "%9.6f" % (f),
            print

        print "MOX Lattice"
        for m in xrange(mox_lat.num_mixtures()):
            vf = mox_lat.f(m)
            print "%4i" % (m),
            for f in vf:
                print "%9.6f" % (f),
            print

    # make the mixtable for the combined lattices by appending the mox table
    # to the UO2 table (don't include the clean mixtures at the front of the 
    # table)
    num_mixtures = uo2_lat.num_mixtures() + mox_lat.num_mixtures() - 4
    table        = Vec_Dbl(num_mixtures * 4)
    ctr          = 0
    mox_offset   = uo2_lat.num_mixtures()

    # add UO2 mixtures
    for m in xrange(uo2_lat.num_mixtures()):
        vf = uo2_lat.f(m)
        for f in vf:
            table[ctr] = f
            ctr        = ctr + 1

    # add MOX mixtures, skipping the clean mixes
    for m in xrange(4, mox_lat.num_mixtures()):
        vf = mox_lat.f(m)
        for f in vf:
            table[ctr] = f
            ctr        = ctr + 1

    # make the cleanids
    cleanids = [0, 1, 2, 3]

    # the total core is 3x3 assemblies (2x2 fuel surrounded by water)
    xylat  = uo2_lat.xy_planes()
    Nr     = len(xylat) - 1
    delta  = Vec_Dbl(Nr, 0.0)
    for i in xrange(Nr):
        delta[i] = xylat[i+1] - xylat[i]

    if Nr % 2 != 0: 
        print "Non-even lattices cells."
        sys.exit(1)

    # build the core planes
    xycore = Vec_Dbl(int(2.5*Nr) + 1, 0.0) 
    for n in xrange(2):
        for i in xrange(Nr):
            index             = i + n * Nr
            xycore[index + 1] = xycore[index] + delta[i]
    for i in xrange(Nr/2):
        index             = i + 2 * Nr
        xycore[index + 1] = xycore[index] + delta[i]

    # z-planes (14 in each assembly)
    height = 14.28 * 1.5
    Nz     = 21
    z      = [0.0] * (Nz + 1)
    dz     = height / float(Nz)
    for k in xrange(Nz):
        z[k+1] = z[k] + dz

    # get matids for each lattice
    uo2ids = Vec_Int(uo2_lat.mixids())
    moxids = Vec_Int(mox_lat.mixids())
    
    # update the mox mixtures (leave clean zones alone)
    for m in xrange(len(moxids)):
        if moxids[m] > 3:
            moxids[m] = moxids[m] + mox_offset - 4

    # assign the matids
    Nx = len(xycore) - 1
    Ny = len(xycore) - 1

    # arrangement
    # |-----|-----|-----|
    # |     |     |     |
    # | mod | mod | mod |
    # |     |     |     |
    # |-----|-----|-----|
    # |     |     |     |
    # | mox | uo2 | mod | y
    # |     |     |     |
    # |-----|-----|-----|
    # |     |     |     |
    # | uo2 | mox | mod |
    # |     |     |     |
    # |-----|-----|-----|
    #          x

    mixids = Vec_Int(Nx * Ny * Nz, 3)

    kend = Nz / 2

    # (0, 0) lattice
    for k in xrange(kend):
        for j in xrange(Nr):
            for i in xrange(Nr):
                lat_cell     = i + j * Nr
                cell         = i + j * Ny + k * Nx * Ny
                mixids[cell] = uo2ids[lat_cell]
    
    # (1, 0) lattice
    for k in xrange(kend):
        for j in xrange(Nr):
            for i in xrange(Nr):
                lat_cell     = i + j * Nr
                cell         = (i + Nr) + j * Ny + k * Nx * Ny
                mixids[cell] = moxids[lat_cell]
    
    # (0, 1) lattice
    for k in xrange(kend):
        for j in xrange(Nr):
            for i in xrange(Nr):
                lat_cell     = i + j * Nr
                cell         = i + (j + Nr) * Ny + k * Nx * Ny
                mixids[cell] = moxids[lat_cell]
    
    # (1, 1) lattice
    for k in xrange(kend):
        for j in xrange(Nr):
            for i in xrange(Nr):
                lat_cell     = i + j * Nr
                cell         = (i + Nr) + (j + Nr) * Ny + k * Nx * Ny
                mixids[cell] = uo2ids[lat_cell]

    return (xycore, z, mixids, cleanids, table)

##---------------------------------------------------------------------------##
## DB
##---------------------------------------------------------------------------##

entries = {
    "problem_type"           : "FIXED_SOURCE",
    "num_groups"             : 7,
    "downscatter"            : False,
    "Pn_order"               : 0,
    "tolerance"              : 1.0e-3,
    "max_itr"                : 400,
    "linear_solver_xml_file" : "azilut04.xml",
    "boundary"               : "reflect",
    "boundary_db"            : {"reflect" : [1, 0, 1, 0, 1, 0]},
    "SPN_order"              : 1
    }

db = DB.from_dict(entries)

# decomposition
if nodes() == 1:

    db.insert("num_blocks_i", 1)
    db.insert("num_blocks_j", 1)

elif nodes() == 2:

    db.insert("num_blocks_i", 2)
    db.insert("num_blocks_j", 1)

elif nodes() == 16:

    db.insert("num_blocks_i", 4)
    db.insert("num_blocks_j", 4)

# Mesh
(r, z, mixids, cleanids, table) = build_mesh(10)
db.insert("x_edges", r)
db.insert("y_edges", r)
db.insert("z_edges", z)

##---------------------------------------------------------------------------##
## MANAGER
##---------------------------------------------------------------------------##

# make manager, material, and angles
manager = Manager()
mat     = Mat()

# partition the problem
manager.partition(db, mat)

# get mapping and mesh objects
mapp    = manager.get_map()
indexer = manager.get_indexer()
mesh    = manager.get_mesh()

# global and local cell numbers
Gx = indexer.num_global(X)
Gy = indexer.num_global(Y)
Gz = mesh.num_cells_dim(Z)
Nx = mesh.num_cells_dim(X)
Ny = mesh.num_cells_dim(Y)
Nz = mesh.num_cells_dim(Z)

if node() == 0:
    print ">>> Partitioned global mesh with %i x %i x %i cells" \
          % (Gx, Gy, Gz)

##---------------------------------------------------------------------------##
## MATERIAL SETUP
##---------------------------------------------------------------------------##

# vacuum    = 0
# UO2       = 1
# MOX       = 2
# moderator = 3

# set database
xsdb = XS_DB(db)
xsdb.set_num(4)

xsdb.assign_zero(0)
for g in xrange(0, xsdb.num_groups()):
    xsdb.assign_upscatter(1, g, T_UO2[g], U_UO2[g], S_UO2[g])
    xsdb.assign_upscatter(2, g, T_MOX43[g], U_MOX43[g], S_MOX43[g])
    xsdb.assign_upscatter(3, g, T_MOD1[g], U_MOD1[g], S_MOD1[g])

## Assign fission data
xsdb.assign_fission(1, NUF_UO2, C_UO2)
xsdb.assign_fission(2, NUF_MOX43, C_MOX43)

# make macro mixer
mixer = Macro_Mixer(xsdb)
mixer.set(cleanids, table)

# make the material database
mixer.mix_with_global_ids(mixids, mat)

##---------------------------------------------------------------------------##
## ENERGY PARTITIONING
##---------------------------------------------------------------------------##

manager.partition_energy(mat)

##---------------------------------------------------------------------------##
## SOURCE SETUP
##---------------------------------------------------------------------------##

# allocate source and problem state
source  = Isotropic_Source()
manager.setup(source)

total = Gx * Gy * Gz
Ng    = mat.num_groups()

srcids     = Vec_Int(total, 0)
srcstr     = Vec_Dbl(total, 0.0)
num_shapes = 2
shapes     = Vec_Dbl(2 * mat.num_groups(), 0.0)

chi0 = xsdb.fission_data(1, 0, CHI)
chi1 = xsdb.fission_data(2, 0, CHI)

# source 0 spectrum -> UO2 Chi
# source 1 spectrum -> MOX Chi

# make shapes
ctr = 0
for g in xrange(Ng):
    shapes[ctr] = xsdb.fission_data(1, g, CHI)
    ctr        += 1
for g in xrange(Ng):
    shapes[ctr] = xsdb.fission_data(2, g, CHI)
    ctr        += 1

# assign ids and strengths
for cell in xrange(total):
    matid = mixids[cell]
    if mat.assigned_fission(matid):
        for g in xrange(Ng):
            srcstr[cell] += mat.fission_data(matid, g, NU_SIGMA_F)
        if mat.fission_data(matid, 0, CHI) == chi1:
            srcids[cell] = 1

# set the source
source.set(num_shapes, shapes, srcids, srcstr)

##---------------------------------------------------------------------------##
## SOLVE
##---------------------------------------------------------------------------##

if node() == 0:
    print ">>> Setup complete"
    print ">>> Solving with %s differencing" % (manager.spatial_descriptor())

# solve the problem
manager.solve(source)

##---------------------------------------------------------------------------##
## OUTPUT
##---------------------------------------------------------------------------##

# make SILO output
silo = SILO()
silo.add_mixer(mixer)
silo.open("fs")

phi = Vec_Dbl(mesh.num_cells(), 0.0)

for g in xrange(Ng):
    flux = manager.moments(g)
    for cell in xrange(mesh.num_cells()):
        phi[cell] = phi[cell] + flux.scalar_flux(cell)
silo.add("phi", phi)

silo.close()

##---------------------------------------------------------------------------##
## TIMING
##---------------------------------------------------------------------------##

# output final database (has class-dependent defaults)
db.output()

timer.stop()
time = timer.wall_clock()

keys = timer_keys()
if len(keys) > 0 and node() == 0:
    print "\n"
    print "TIMING : Problem ran in %16.6e seconds." % (time)
    print "------------------------------------------------------------------"
    for key in keys:
        print "%30s : %16.6e %16.6e" % (key, timer_value(key) / time, timer_value(key))
    print "------------------------------------------------------------------"

##---------------------------------------------------------------------------##

manager.close()
finalize()

###############################################################################
## end of fs.py
###############################################################################
