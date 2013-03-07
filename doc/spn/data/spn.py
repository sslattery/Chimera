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
num_groups = 10
upscatter = False
##---------------------------------------------------------------------------##


entries = {
    "delta_x"     : 1.0,
    "delta_y"     : 1.0,
    "delta_z"     : 1.0,
    "num_cells_i" : 5,
    "num_cells_j" : 5,
    "num_cells_k" : 5,
    "num_groups"  : num_groups,
    "SPN_order"   : spn_order,
    "boundary"    : "reflect",
    "boundary_db" : {"reflect" : [1,1,1,1,1,1]},
    "Pn_order"    : pn_order,
    "downscatter" : True
    }

db = DB.from_dict(entries)

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
    ingroup_xs.append(0.25)
    down_xs.append(1.0)

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

dg0 = [ingroup_xs]
dg1 = [down_xs,    ingroup_xs]
dg2 = [down_xs,    down_xs,    ingroup_xs]
dg3 = [down_xs,    down_xs,    down_xs,    ingroup_xs]
dg4 = [down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]
dg5 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]
dg6 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]
dg7 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]
dg8 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]
dg9 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]

mdown = [dg0, dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9]
print dsxs
print mdown
# upscatter cross sections.
up_xs = []
for m in xrange(pn_order+1):
    up_xs.append(0.1)

usxs = []
cup = []
local_cup = []
for g in xrange(num_groups):
    for j in xrange(num_groups-g-1):
        local_cup.append(j)
    cup.append(local_cup)
    local_cup = []
    for i in xrange(num_groups):
        if i == g:
            usxs.append(ingroup_xs)
        elif i < g :
            usxs.append(down_xs)
        else:
            usxs.append(up_xs)

ug0 = [ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs]
ug1 = [down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs]
ug2 = [down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs]
ug3 = [down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs,      up_xs,      up_xs]
ug4 = [down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs,      up_xs]
ug5 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs,      up_xs]
ug6 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs,      up_xs]
ug7 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs,      up_xs]
ug8 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs, up_xs]
ug9 = [down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    down_xs,    ingroup_xs]

mup = [ug0, ug1, ug2, ug3, ug4, ug5, ug6, ug7, ug8, ug9]
cup = [[1,2,3,4,5,6,7,8,9],[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7],[1,2,3,4,5,6],[1,2,3,4,5],[1,2,3,4],[1,2,3],[1,2],[1],[]]

# assign xs
sigma_t = 2.0

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

##---------------------------------------------------------------------------##

manager.close()
finalize()
