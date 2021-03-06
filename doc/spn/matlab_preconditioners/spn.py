import sys
import math

from spn_fv import *

##---------------------------------------------------------------------------##

initialize(sys.argv)


##---------------------------------------------------------------------------##
# SPECTRAL ANALYSIS PARAMETERS
##---------------------------------------------------------------------------##
pn_order = 1
spn_order = 1
num_groups = 1
upscatter = True

sigma_t = 5.0
sigma_ing = 2.0
sigma_down = 3.0
sigma_up = 1.0

##---------------------------------------------------------------------------##


entries = {
    "delta_x"     : 0.03,
    "delta_y"     : 0.03,
    "delta_z"     : 0.03,
    "num_cells_i" : 6,
    "num_cells_j" : 6,
    "num_cells_k" : 6,
    "num_groups"  : num_groups,
    "SPN_order"   : spn_order,
    "boundary"    : "reflect",
    "boundary_db" : {"reflect" : [1,1,1,1,1,1]},
    "Pn_order"    : pn_order,
    "downscatter" : False
    }

db = DB.from_dict(entries)

db.insert("linear_solver_xml_file","mcls.xml")

## DECOMPOSITION

if nodes() == 1:

    db.insert("num_blocks_i", 1)
    db.insert("num_blocks_j", 1)

elif nodes() == 2:

    db.insert("num_blocks_i", 2)
    db.insert("num_blocks_j", 1)

elif nodes() == 4:

    db.insert("num_blocks_i", 2)
    db.insert("num_blocks_j", 2)

##---------------------------------------------------------------------------##

manager = Manager()
mat     = Mat()
source  = General_Source()

manager.partition(db, mat)

mesh = manager.get_mesh()

##---------------------------------------------------------------------------##

mat.set_num(1)
mat.assign_id(0, [])

# downscatter cross sections
ingroup_xs = []
down_xs = []
for m in xrange(pn_order+1):
    ingroup_xs.append(sigma_ing)
    down_xs.append(sigma_down)

dsxs = []
local_dsxs = []
cdown = []
for g in xrange(num_groups):
    cdown.append([])
    for i in xrange(g+1):
        if i == g:
            local_dsxs.append(ingroup_xs)
        else:
            local_dsxs.append(down_xs)
    dsxs.append(local_dsxs)
    local_dsxs = []

# upscatter cross sections
up_xs = []
for m in xrange(pn_order+1):
    up_xs.append(sigma_up)

usxs = []
local_usxs = []
cup = []
local_cup = []
for g in xrange(num_groups):
    for j in xrange(num_groups-g-1):
        local_cup.append(j+g+1)
    cup.append(local_cup)
    local_cup = []
    for i in xrange(num_groups):
        if i == g:
            local_usxs.append(ingroup_xs)
        elif i < g :
            local_usxs.append(down_xs)
        else:
            local_usxs.append(up_xs)
    usxs.append(local_usxs)
    local_usxs = []

# assign xs
if num_groups == 1:
    mat.assign_xs(0, 0, sigma_t, dsxs[0])
elif upscatter:
    for g in xrange(num_groups):
        mat.assign_upscatter(0, g, sigma_t, cup[g], usxs[g])
else:
    for g in xrange(num_groups):
        mat.assign_upscatter(0, g, sigma_t, cdown[g], dsxs[g])

##---------------------------------------------------------------------------##

manager.partition_energy(mat)

##---------------------------------------------------------------------------##
## External source setup

manager.setup(source)

source.set_num(1)

source.assign_id(0, [])
for g in xrange(num_groups):
    source.assign_isotropic(0, g, 1.0)

##---------------------------------------------------------------------------##

manager.write_Matlab()
#manager.solve(source)

##---------------------------------------------------------------------------##

adapter = manager.get_flux_adapter()

flux = Vec_Dbl()
adapter.get_local_scalar_flux(0, flux)

for i in xrange( mesh.num_cells() ):
    print flux[i]

##---------------------------------------------------------------------------##

manager.close()
finalize()
