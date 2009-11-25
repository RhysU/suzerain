#!/usr/bin/env sage
# Symbolic computations for the test cases in test_timestepper.cpp
#
# Examining particular results can be done by starting sage
# from this directory, typing 'attach test_timestepper.sage'.

# SMR91 coefficients
gamma = [8/15, 5/12, 3/4];
zeta  = [0, -17/60, -5/12];
alpha = [29/96, -3/40, 1/6];
beta  = [37/160, 5/24, 1/6];

# Purely explicit Riccati equation nonlinear operator
# is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
def RiccatiNonlinearOperator(a, b):
    return lambda state: map(lambda x: x**2 + b*x - a**2 - a*b ,state)

def substep_explicit(op, a, b, delta_t, substep):
    for i in range(0,len(a)):
        b[i] = delta_t * zeta[substep] * b[i] + a[i]
        a[i]= op([a[i]])[0]
        b[i] = b[i] + delta_t * gamma[substep] * a[i]
    return (a, b)

substep_explicit_0 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 0)
substep_explicit_1 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 1)
substep_explicit_2 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 2)
