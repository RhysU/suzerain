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
eta   = [0, alpha[0] + beta[0], alpha[0] + beta[0] + alpha[1] + beta[1]];

# Purely explicit Riccati equation nonlinear operator
# is the right hand side of (d/dt) y = y^2 + b y - a^2 -a b
def RiccatiNonlinearOperator(a, b):
    return lambda y, t: y**2 + b*y - a**2 - a*b

def substep_explicit(op, t, a, b, delta_t, substep):
    for i in range(0,len(a)):
        b[i] = delta_t * zeta[substep] * b[i] + a[i]
        a[i]= op(a[i], t + eta[substep] * delta_t)
        b[i] = b[i] + delta_t * gamma[substep] * a[i]
    return (a, b)

# Manufactured answers for the time-independent, explicit substep test case
substep_explicit_0 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), 0, [5, 7], [11, 13], 17, 0)
substep_explicit_1 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), 0, [5, 7], [11, 13], 17, 1)
substep_explicit_2 = substep_explicit(
    RiccatiNonlinearOperator(2, 3), 0, [5, 7], [11, 13], 17, 2)

# Purely explicit operator for (d/dt) = cos(t)
def CosineNonlinearOperator():
    return lambda y, t: cos(t)

# Manufactured answers for the time-independent, explicit substep test case
substep_explicit_0_td = substep_explicit(
    CosineNonlinearOperator(), pi/3, [5, 7], [11, 13], 17, 0)
substep_explicit_1_td = substep_explicit(
    CosineNonlinearOperator(), pi/3, [5, 7], [11, 13], 17, 1)
substep_explicit_2_td = substep_explicit(
    CosineNonlinearOperator(), pi/3, [5, 7], [11, 13], 17, 2)

# Nonlinear portion of the Riccati operator is the right hand side of (d/dt) y
# = y^2 + b y - a^2 -a b minus the b y term
def RiccatiHybridNonlinearOperator(a, b):
    return lambda y: y**2 - a**2 - a*b

# Linear portion of the Riccati operator
def RiccatiHybridLinearOperator(a, b):
    return lambda y: b*y

def substep_hybrid(linear, nonlinear, a, b, delta_t, substep):
    linear_coefficient = linear(1)
    for i in range(0,len(a)):
        b[i] = (delta_t * zeta[substep] * b[i]
                    + a[i] + delta_t*alpha[substep]*linear(a[i]))
        a[i] = nonlinear(a[i])
        b[i] = b[i] + delta_t * gamma[substep] * a[i]
        b[i] = (1 - delta_t * beta[substep] * linear_coefficient)^(-1)*b[i]
    return (a, b)

# Manufactured answers for the hybrid substep test case
substep_hybrid_0 = substep_hybrid(
    RiccatiHybridLinearOperator(2, 3),
    RiccatiHybridNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 0)
substep_hybrid_1 = substep_hybrid(
    RiccatiHybridLinearOperator(2, 3),
    RiccatiHybridNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 1)
substep_hybrid_2 = substep_hybrid(
    RiccatiHybridLinearOperator(2, 3),
    RiccatiHybridNonlinearOperator(2, 3), [5, 7], [11, 13], 17, 2)

# Checking the residual for the hybrid substep test case
def substep_hybrid_residual(linear, nonlinear,
                            a_old, b_old, a_new, b_new,
                            delta_t, substep):
    retval = []
    for i in range(0,len(a_old)):
        lhs = (
                b_new[i]
                - delta_t * beta[substep] * linear(b_new[i])
            )
        rhs = (
                a_old[i]
                + delta_t * alpha[substep] * linear(a_old[i])
                + delta_t * gamma[substep] * a_new[i]
                + delta_t * zeta[substep]  * b_old[i]
            )
        retval.append(lhs - rhs)
    return retval

# Ensure the manufactured answers solve the original equation
substep_hybrid_residual_0 = substep_hybrid_residual(
        RiccatiHybridLinearOperator(2,3),
        RiccatiHybridNonlinearOperator(2,3),
        [5,7], [11,13],
        substep_hybrid_0[0], substep_hybrid_0[1],
        17, 0
    )
substep_hybrid_residual_1 = substep_hybrid_residual(
        RiccatiHybridLinearOperator(2,3),
        RiccatiHybridNonlinearOperator(2,3),
        [5,7], [11,13],
        substep_hybrid_1[0], substep_hybrid_1[1],
        17, 1
    )
substep_hybrid_residual_2 = substep_hybrid_residual(
        RiccatiHybridLinearOperator(2,3),
        RiccatiHybridNonlinearOperator(2,3),
        [5,7], [11,13],
        substep_hybrid_2[0], substep_hybrid_2[1],
        17, 2
    )
