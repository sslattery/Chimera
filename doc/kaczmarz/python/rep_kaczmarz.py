import numpy
import math
import random

m = 3

# Kaczmarz projection method
def kaczmarz(A,x,b):
    r = b - A*x
    num_iters = 0
    while numpy.linalg.norm(r,numpy.inf)/numpy.linalg.norm(b,numpy.inf) > 1.0e-2:
        i = random.randint(0,m)
        j = m-i-1
        alpha = numpy.asscalar( numpy.dot(A[j],A[j].T) )
        scalarb = numpy.asscalar(b[j])
        scalarAjdotx = numpy.asscalar( numpy.dot(A[j],x) )
        x = x + (1/alpha)*(scalarb - scalarAjdotx)*(A.T)[j].T
        r = b - A*x
        num_iters += 1
    print "Kaczmarz iterations to converge:", num_iters
    return x

# Solve a 3x3 linear system Ax=b with initial state x_i = 0
A = numpy.matrix([ [1.0, -0.3, -0.4],[-0.2, 1.0, -0.4],[-0.5, -0.4, 1.0] ])
b = numpy.matrix([ [10.0], [1.0], [1.0] ])
x_kaczmarz_1 = numpy.matrix([ [0.0], [0.0], [0.0] ])
x_kaczmarz_2 = numpy.matrix([ [0.0], [0.0], [0.0] ])
x_kaczmarz_3 = numpy.matrix([ [0.0], [0.0], [0.0] ])
x_kaczmarz_4 = numpy.matrix([ [0.0], [0.0], [0.0] ])
x_adjoint = numpy.matrix([ [0.0], [0.0], [0.0] ])

# Spectral Radius of A
eigenvalues_A, eigenvectors_A = numpy.linalg.eig(A)
rho_A = max(eigenvalues_A)
print 
print "Spectral radius A =", rho_A
print

# Kaczmarz solution 1
x_kaczmarz_1 = kaczmarz(A,x_kaczmarz_1,b)
print
print "Kaczmarz"
print x_kaczmarz_1

# Kaczmarz solution 2
x_kaczmarz_2 = kaczmarz(A,x_kaczmarz_2,b)
print
print "Kaczmarz"
print x_kaczmarz_2

# Kaczmarz solution 3
x_kaczmarz_3 = kaczmarz(A,x_kaczmarz_3,b)
print
print "Kaczmarz"
print x_kaczmarz_3

# Kaczmarz solution 4
x_kaczmarz_4 = kaczmarz(A,x_kaczmarz_4,b)
print
print "Kaczmarz"
print x_kaczmarz_4

# Reference solution
x_numpy = numpy.linalg.solve( A, b )
print
print "Reference"
print x_numpy

