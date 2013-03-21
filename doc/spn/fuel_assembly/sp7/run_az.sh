#!/bin

#
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut01.py > azilut01.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut02.py > azilut02.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut03.py > azilut03.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut04.py > azilut04.log
