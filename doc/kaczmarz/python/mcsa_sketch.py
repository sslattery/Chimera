import numpy
import math
import random
import time

##---------------------------------------------------------------------------##
# compute the inf norm of a vector
def inf_norm(x):
    m = len(x)
    error = []
    for i in xrange(m):
        error.append( x[i][0] )

    return max(error)

##---------------------------------------------------------------------------##
# Sampling function for stochastic matrices ( O(n) algorithm )
def sample(P,s):
    cdf = 0.0
    zeta = random.random()
    inverted = 0
    for i in xrange( len(P[s]) ):
        cdf += P[s][i]
        if zeta < cdf and not inverted:
            ns = i
            inverted = 1

    return int(ns)

##---------------------------------------------------------------------------##
# Direct Monte Carlo solution method
def direct_mc(np,wc,H,P,b,x):
    m = len(x)
    for i in xrange(m):
        for n in xrange(np):
            s = i
            walk = 1
            w = 1.0
            while walk:
                x[i][0] += w*b[s][0]
                ns = sample(P,s)
                if P[s][ns] == 0:
                    w = 0
                else:
                    w *= H[s][ns]/P[s][ns]
                if w < wc:
                    walk = 0
                s = ns
            
        x[i][0] /= np

##---------------------------------------------------------------------------##
# Adjoint Monte Carlo solution method.
def adjoint_mc(np,wc,H,Q,b,x):
    m = len(x)

    b_norm = 0.0
    for j in xrange(m):
        b_norm += b[j][0]

    for n in xrange(np):
        init = 0
        b_cdf = 0
        zeta = random.random()
        inverted = 0
        for j in xrange(m):
            b_cdf += b[j][0]/b_norm
            if zeta < b_cdf and not inverted:
                init = j
                inverted = 1
        walk = 1
        w = b_norm / b[init][0]
        wf = wc*w
        s = init
        while walk:
            x[s][0] += w*b[init][0]
            ns = sample(Q,s)
            if Q[s][ns] == 0:
                w = 0
            else:
                w *= H[ns][s]/Q[s][ns]
            if w < wf:
                walk = 0
            s = ns

    for i in xrange(m):
        x[i][0] /= np

##---------------------------------------------------------------------------##
# Sequential Monte Carlo solution method with direct residual solve.
def seq_direct_mc(np,wc,A,H,P,b,x,tol):
    m = len(x)
    r = []
    delta_x = []
    for i in xrange(m):
        r.append( [0.0] )
        delta_x.append( [1.0] )
    iterations = 0
    while inf_norm(delta_x) > tol:
        iterations += 1
        for i in xrange(m):
            delta_x[i][0] = 0.0
            Ax_row = 0.0
            for j in xrange(m):
                Ax_row += A[i][j]*x[j][0]
            r[i][0] = b[i][0] - Ax_row
        direct_mc(np,wc,H,P,r,delta_x)
        for i in xrange(m):
            x[i][0] += delta_x[i][0]

    print "Sequential direct iterations: ", iterations

##---------------------------------------------------------------------------##
# Sequential Monte Carlo solution method with adjoint residual solve.
def seq_adjoint_mc(np,wc,A,H,Q,b,x,tol):
    m = len(x)
    r = []
    delta_x = []
    for i in xrange(m):
        r.append( [0.0] )
        delta_x.append( [1.0] )
    iterations = 0
    while inf_norm(delta_x) > tol:
        iterations += 1
        for i in xrange(m):
            delta_x[i][0] = 0.0
            Ax_row = 0.0
            for j in xrange(m):
                Ax_row += A[i][j]*x[j][0]
            r[i][0] = b[i][0] - Ax_row
        adjoint_mc(np,wc,H,Q,r,delta_x)
        for i in xrange(m):
            x[i][0] += delta_x[i][0]

    print "Sequential adjoint iterations:" , iterations

##---------------------------------------------------------------------------##
# Monte Carlo synthetic acceleration with direct solve.
def mcsa_direct(np,wc,A,H,P,b,x,tol):
    m = len(x)
    r = []
    delta_x = []
    for i in xrange(m):
        r.append( [0.0] )
        delta_x.append( [1.0] )
    iterations = 0
    while inf_norm(delta_x) > tol:
        iterations += 1
        for i in xrange(m):
            delta_x[i][0] = 0.0
            ImAx_row = 0.0
            for j in xrange(m):
                if ( i == j ):
                    ImAx_row += (1.0-A[i][j])*x[j][0]
                else:
                    ImAx_row += -A[i][j]*x[j][0]
            x[i][0] += ImAx_row + b[i][0]
        for i in xrange(m):
            Ax_row = 0.0
            for j in xrange(m):
                Ax_row += A[i][j]*x[j][0]
            r[i][0] = b[i][0] - Ax_row
        direct_mc(np,wc,H,P,r,delta_x)
        for i in xrange(m):
            x[i][0] += delta_x[i][0]

    print "MCSA direct iterations:       ", iterations

##---------------------------------------------------------------------------##
# Monte Carlo synthetic acceleration with adjoint solve.
def mcsa_adjoint(np,wc,A,H,Q,b,x,tol):
    m = len(x)
    r = []
    delta_x = []
    for i in xrange(m):
        r.append( [0.0] )
        delta_x.append( [1.0] )
    iterations = 0
    while inf_norm(delta_x) > tol:
        iterations += 1
        for i in xrange(m):
            delta_x[i][0] = 0.0
            ImAx_row = 0.0
            for j in xrange(m):
                if ( i == j ):
                    ImAx_row += (1.0-A[i][j])*x[j][0]
                else:
                    ImAx_row += -A[i][j]*x[j][0]
            x[i][0] += ImAx_row + b[i][0]
        for i in xrange(m):
            Ax_row = 0.0
            for j in xrange(m):
                Ax_row += A[i][j]*x[j][0]
            r[i][0] = b[i][0] - Ax_row
        adjoint_mc(np,wc,H,Q,r,delta_x)
        for i in xrange(m):
            x[i][0] += delta_x[i][0]

    print "MCSA adjoint iterations:      ", iterations

##---------------------------------------------------------------------------##
# Solve a 3x3 linear system Ax=b with initial state x_i = 0
m = 3
A = [ [1.0, -0.3, -0.4],[-0.2, 1.0, -0.4],[-0.5, -0.4, 1.0] ]
b = [ [10.0], [1.0], [1.0] ]
x_direct = [ [0.0], [0.0], [0.0] ]
x_adjoint = [ [0.0], [0.0], [0.0] ]
x_seq_direct = [ [0.0], [0.0], [0.0] ]
x_seq_adjoint = [ [0.0], [0.0], [0.0] ]
x_mcsa_direct = [ [0.0], [0.0], [0.0] ]
x_mcsa_adjoint = [ [0.0], [0.0], [0.0] ]

# Identity matrix
I = [ [1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0] ]

# Iteration matrix
H = []
for i in xrange(m):
    H_row = []
    for j in xrange(m):
        H_row.append( I[i][j]-A[i][j] )
    H.append(H_row)

# Spectral Radius of H
eigenvalues_H, eigenvectors_H = numpy.linalg.eig(H)
rho_H = max(eigenvalues_H)
print 
print "Spectral radius H =", rho_H
print

# Probability Matrix
P = []
for i in xrange(m):
    P_row = []
    P_row_sum = 0.0
    for j in xrange(m):
        P_row_sum += abs(H[i][j])
    for j in xrange(m):
        P_row.append( abs(H[i][j]) / P_row_sum )
    P.append(P_row)

# Adjoint Probability Matrix
Q = []
for i in xrange(m):
    Q_row = []
    Q_col_sum = 0.0
    for j in xrange(m):
        Q_col_sum += abs(H[j][i])
    for j in xrange(m):
        Q_row.append( abs(H[j][i]) / Q_col_sum )
    Q.append(Q_row)

# Setup
np = 100
wc = 1.0E-12
tol = 1.0E-8
print "Histories             =", np
print "Weight cutoff         =", wc
print "Convergence tolerance =", tol
print

# Direct Monte Carlo Solution
begin = time.clock()
direct_mc(np, wc, H, P, b, x_direct)
end = time.clock()
print "Direct time:", end-begin

# Adjoint Monte Carlo Solution
begin = time.clock()
adjoint_mc(np, wc, H, Q, b, x_adjoint)
end = time.clock()
print "Adjoint time:", end-begin

# Sequential Monte Carlo with Direct Solve
begin = time.clock()
seq_direct_mc(np, wc, A, H, P, b, x_seq_direct, tol)
end = time.clock()
print "Sequential Direct time:", end-begin

# Sequential Monte Carlo with Adjoint Solve
begin = time.clock()
seq_adjoint_mc(np, wc, A, H, Q, b, x_seq_adjoint, tol)
end = time.clock()
print "Sequential Adjoint time:", end-begin

# MCSA with Direct Solve
begin = time.clock()
mcsa_direct(np, wc, A, H, P, b, x_mcsa_direct, tol)
end = time.clock()
print "MCSA Direct time:", end-begin

# MCSA with Adjoint Solve
begin = time.clock()
mcsa_adjoint(np, wc, A, H, Q, b, x_mcsa_adjoint, tol)
end = time.clock()
print "MCSA Adjoint time:", end-begin

# Numpy linear solution
A_matrix = numpy.mat(A)
b_matrix = numpy.mat(b)
x_numpy = numpy.linalg.solve( A_matrix, b_matrix )

# Compare the solutions with the inf norm and L2 norm of the relative error
direct_error = []
direct_l2_norm = 0.0
adjoint_error = []
adjoint_l2_norm = 0.0
seq_direct_error = []
seq_direct_l2_norm = 0.0
seq_adjoint_error = []
seq_adjoint_l2_norm = 0.0
mcsa_direct_error = []
mcsa_direct_l2_norm = 0.0
mcsa_adjoint_error = []
mcsa_adjoint_l2_norm = 0.0
for i in xrange(m):
    direct_error.append( \
        abs(x_numpy[i][0]-x_direct[i][0])/abs(x_direct[i][0]) )
    direct_l2_norm += direct_error[i]**2

    adjoint_error.append( \
        abs(x_numpy[i][0]-x_adjoint[i][0])/abs(x_adjoint[i][0]) )
    adjoint_l2_norm += adjoint_error[i]**2

    seq_direct_error.append( \
        abs(x_numpy[i][0]-x_seq_direct[i][0])/abs(x_seq_direct[i][0]) )
    seq_direct_l2_norm += seq_direct_error[i]**2

    seq_adjoint_error.append( \
        abs(x_numpy[i][0]-x_seq_adjoint[i][0])/abs(x_seq_adjoint[i][0]) )
    seq_adjoint_l2_norm += seq_adjoint_error[i]**2

    mcsa_direct_error.append( \
        abs(x_numpy[i][0]-x_mcsa_direct[i][0])/abs(x_mcsa_direct[i][0]) )
    mcsa_direct_l2_norm += mcsa_direct_error[i]**2

    mcsa_adjoint_error.append( \
        abs(x_numpy[i][0]-x_mcsa_adjoint[i][0])/abs(x_mcsa_adjoint[i][0]) )
    mcsa_adjoint_l2_norm += mcsa_adjoint_error[i]**2

direct_inf_norm = max(direct_error)
direct_l2_norm = math.sqrt(direct_l2_norm)

adjoint_inf_norm = max(adjoint_error)
adjoint_l2_norm = math.sqrt(adjoint_l2_norm)

seq_direct_inf_norm = max(seq_direct_error)
seq_direct_l2_norm = math.sqrt(seq_direct_l2_norm)

seq_adjoint_inf_norm = max(seq_adjoint_error)
seq_adjoint_l2_norm = math.sqrt(seq_adjoint_l2_norm)

mcsa_direct_inf_norm = max(mcsa_direct_error)
mcsa_direct_l2_norm = math.sqrt(mcsa_direct_l2_norm)

mcsa_adjoint_inf_norm = max(mcsa_adjoint_error)
mcsa_adjoint_l2_norm = math.sqrt(mcsa_adjoint_l2_norm)

# Output
print
print "Direct MC Solution"
print x_direct
print

print "Adjoint MC Solution"
print x_adjoint
print

print "Sequential Direct MC Solution"
print x_seq_direct
print

print "Sequential Adjoint MC Solution"
print x_seq_adjoint
print

print "MCSA Direct Solution"
print x_mcsa_direct
print

print "MCSA Adjoint Solution"
print x_mcsa_adjoint
print

print "Numpy Solution"
print x_numpy
print

print "Direct inf_norm             =", direct_inf_norm
print "Direct L2_norm              =", direct_l2_norm
print

print "Adjoint inf_norm            =", adjoint_inf_norm
print "Adjoint L2_norm             =", adjoint_l2_norm
print

print "Sequential Direct inf_norm  =", seq_direct_inf_norm
print "Sequential Direct L2_norm   =", seq_direct_l2_norm
print

print "Sequential Adjoint inf_norm =", seq_adjoint_inf_norm
print "Sequential Adjoint L2_norm  =", seq_adjoint_l2_norm
print

print "MCSA Direct inf_norm        =", mcsa_direct_inf_norm
print "MCSA Direct L2_norm         =", mcsa_direct_l2_norm
print

print "MCSA Adjoint inf_norm       =", mcsa_adjoint_inf_norm
print "MCSA Adjoint L2_norm        =", mcsa_adjoint_l2_norm
print
