#!/bin

mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilu00.py > azilu00.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilu01.py > azilu01.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilu02.py > azilu02.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilu10.py > azilu10.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilu20.py > azilu20.log
#
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut01.py > azilut01.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut02.py > azilut02.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut03.py > azilut03.log
mpirun -np 16 -machinefile ~/.machines ./pyspn fs_azilut04.py > azilut04.log
