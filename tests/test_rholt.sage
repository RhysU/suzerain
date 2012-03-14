#!/usr/bin/env sage
# Symbolic computations for the test cases in test_rholt.cpp
#
# Examining particular results can be done by starting sage
# from this directory, typing 'attach test_rholt.sage',
# and then examining high precision results using, e.g.
#   e.subs(x=1,y=2,z=3).n(200);
#   u[0].subs(x=1,y=2,z=3).n(200);
#   tau[1][2].subs(x=1,y=2,z=3).n(200);

# Variable and field declarations
x = var('x');
y = var('y');
z = var('z');

# Constants, including non-zero bulk viscosity
gamma = 14/10;
beta  = 2/3;
alpha = 5;

# Conservative state variables
rho = 2*(x^2)*y*z + 3*x*(y^2)*z + 5*x*y*(z^2);
m = [
    7*sin( 7*x) + 11*sin(11*y) + 13*sin(13*z),
    17*sin(17*x) + 19*sin(19*y) + 23*sin(23*z),
    31*sin(31*x) + 37*sin(37*y) + 41*sin(41*z)
];
e = 311*(x^2)*y*z + 313*x*(y^2)*z + 317*x*y*(z^2);

# Gradients of conservative state variables
grad_rho = [ rho.diff(x), rho.diff(y), rho.diff(z) ];
div_grad_rho = grad_rho[0].diff(x) + grad_rho[1].diff(y) + grad_rho[2].diff(z);
grad_grad_rho = [
    [ grad_rho[0].diff(x), grad_rho[0].diff(y), grad_rho[0].diff(z) ],
    [ grad_rho[1].diff(x), grad_rho[1].diff(y), grad_rho[1].diff(z) ],
    [ grad_rho[2].diff(x), grad_rho[2].diff(y), grad_rho[2].diff(z) ]
];
grad_m = [
    [ m[0].diff(x), m[0].diff(y), m[0].diff(z) ],
    [ m[1].diff(x), m[1].diff(y), m[1].diff(z) ],
    [ m[2].diff(x), m[2].diff(y), m[2].diff(z) ]
];
grad_m_grad_m_T = [
    [ grad_m[0][0] + grad_m[0][0],
      grad_m[0][1] + grad_m[1][0],
      grad_m[0][2] + grad_m[2][0] ],
    [ grad_m[1][0] + grad_m[0][1],
      grad_m[1][1] + grad_m[1][1],
      grad_m[1][2] + grad_m[2][1] ],
    [ grad_m[2][0] + grad_m[0][2],
      grad_m[2][1] + grad_m[1][2],
      grad_m[2][2] + grad_m[2][2] ]
];
div_m = m[0].diff(x) + m[1].diff(y) + m[2].diff(z);
grad_div_m = [ div_m.diff(x), div_m.diff(y), div_m.diff(z) ];
div_grad_m = [
    m[0].diff(x,x) + m[0].diff(y,y) + m[0].diff(z,z),
    m[1].diff(x,x) + m[1].diff(y,y) + m[1].diff(z,z),
    m[2].diff(x,x) + m[2].diff(y,y) + m[2].diff(z,z),
];
grad_e = [ e.diff(x), e.diff(y), e.diff(z) ];
div_grad_e = grad_e[0].diff(x) + grad_e[1].diff(y) + grad_e[2].diff(z);

# Classical state variables
p = (gamma - 1)*(e - (m[0]*m[0]+m[1]*m[1]+m[2]*m[2])/(2*rho));
T = gamma*p/rho;
mu = T^beta;
l = (alpha-2/3)*mu;

# Gradients of classical state variables
grad_p  = [  p.diff(x),  p.diff(y),  p.diff(z) ];
grad_T  = [  T.diff(x),  T.diff(y),  T.diff(z) ];
grad_mu = [ mu.diff(x), mu.diff(y), mu.diff(z) ];
grad_l  = [  l.diff(x),  l.diff(y),  l.diff(z) ];

# Second derivatives of classical state variables
div_grad_p  = grad_p[0].diff(x)  + grad_p[1].diff(y)  + grad_p[2].diff(z);
div_grad_T  = grad_T[0].diff(x)  + grad_T[1].diff(y)  + grad_T[2].diff(z);
div_grad_mu = grad_mu[0].diff(x) + grad_mu[1].diff(y) + grad_mu[2].diff(z);
div_grad_l  = grad_l[0].diff(x)  + grad_l[1].diff(y)  + grad_l[2].diff(z);

# Velocity field and its derivatives
u = [
    m[0]/rho,
    m[1]/rho,
    m[2]/rho
];
div_u = u[0].diff(x) + u[1].diff(y) + u[2].diff(z);
grad_u = [
    [ u[0].diff(x), u[0].diff(y), u[0].diff(z) ],
    [ u[1].diff(x), u[1].diff(y), u[1].diff(z) ],
    [ u[2].diff(x), u[2].diff(y), u[2].diff(z) ]
];
div_grad_u = [
    u[0].diff(x,x) + u[0].diff(y,y) + u[0].diff(z,z),
    u[1].diff(x,x) + u[1].diff(y,y) + u[1].diff(z,z),
    u[2].diff(x,x) + u[2].diff(y,y) + u[2].diff(z,z),
];
grad_div_u = [
    div_u.diff(x),
    div_u.diff(y),
    div_u.diff(z)
];
grad_u_grad_u_T = [
    [ grad_u[0][0] + grad_u[0][0],
      grad_u[0][1] + grad_u[1][0],
      grad_u[0][2] + grad_u[2][0] ],
    [ grad_u[1][0] + grad_u[0][1],
      grad_u[1][1] + grad_u[1][1],
      grad_u[1][2] + grad_u[2][1] ],
    [ grad_u[2][0] + grad_u[0][2],
      grad_u[2][1] + grad_u[1][2],
      grad_u[2][2] + grad_u[2][2] ]
];

# Viscous stress tensor
tau = [
    [
        mu*(grad_u[0][0] + grad_u[0][0]) + l*div_u,
        mu*(grad_u[0][1] + grad_u[1][0]),
        mu*(grad_u[0][2] + grad_u[2][0])
    ],
    [
        mu*(grad_u[1][0] + grad_u[0][1]),
        mu*(grad_u[1][1] + grad_u[1][1]) + l*div_u,
        mu*(grad_u[1][2] + grad_u[2][1])
    ],
    [
        mu*(grad_u[2][0] + grad_u[0][2]),
        mu*(grad_u[2][1] + grad_u[1][2]),
        mu*(grad_u[2][2] + grad_u[2][2]) + l*div_u
    ]
];

# Other quantities
div_e_u = (e*u[0]).diff(x) + (e*u[1]).diff(y) + (e*u[2]).diff(z);
div_mu_grad_T = (
      (mu*grad_T[0]).diff(x)
    + (mu*grad_T[1]).diff(y)
    + (mu*grad_T[2]).diff(z)
);
div_p_u = (p*u[0]).diff(x) + (p*u[1]).diff(y) + (p*u[2]).diff(z);
div_tau = [
    tau[0][0].diff(x) + tau[1][0].diff(y) + tau[2][0].diff(z),
    tau[0][1].diff(x) + tau[1][1].diff(y) + tau[2][1].diff(z),
    tau[0][2].diff(x) + tau[1][2].diff(y) + tau[2][2].diff(z),
]
tau_u = [
    tau[0][0]*u[0] + tau[0][1]*u[1] + tau[0][2]*u[2],
    tau[1][0]*u[0] + tau[1][1]*u[1] + tau[1][2]*u[2],
    tau[2][0]*u[0] + tau[2][1]*u[1] + tau[2][2]*u[2],
];
div_tau_u = tau_u[0].diff(x) + tau_u[1].diff(y) + tau_u[2].diff(z);
m_dot_m = m[0]*m[0] + m[1]*m[1] + m[2]*m[2]
grad_m_m = [
    grad_m[0][0]*m[0] + grad_m[0][1]*m[1] + grad_m[0][2]*m[2],
    grad_m[1][0]*m[0] + grad_m[1][1]*m[1] + grad_m[1][2]*m[2],
    grad_m[2][0]*m[0] + grad_m[2][1]*m[1] + grad_m[2][2]*m[2],
];
grad_m_T_m = [
    grad_m[0][0]*m[0] + grad_m[1][0]*m[1] + grad_m[2][0]*m[2],
    grad_m[0][1]*m[0] + grad_m[1][1]*m[1] + grad_m[2][1]*m[2],
    grad_m[0][2]*m[0] + grad_m[1][2]*m[1] + grad_m[2][2]*m[2],
];
u_dot_mu_div_grad_u = (
      u[0]*mu*div_grad_u[0]
    + u[1]*mu*div_grad_u[1]
    + u[2]*mu*div_grad_u[2]
);
u_dot_mu_plus_lambda_grad_div_u = (
      u[0]*(mu+l)*grad_div_u[0]
    + u[1]*(mu+l)*grad_div_u[1]
    + u[2]*(mu+l)*grad_div_u[2]
);

# Divergence form of the convective derivative
u_outer_m = [
    [ u[0]*m[0], u[0]*m[1], u[0]*m[2] ],
    [ u[1]*m[0], u[1]*m[1], u[1]*m[2] ],
    [ u[2]*m[0], u[2]*m[1], u[2]*m[2] ]
];
div_u_outer_m = [
    (
          u_outer_m[0][0].diff(x)
        + u_outer_m[1][0].diff(y)
        + u_outer_m[2][0].diff(z)
    ),
    (
          u_outer_m[0][1].diff(x)
        + u_outer_m[1][1].diff(y)
        + u_outer_m[2][1].diff(z)
    ),
    (
          u_outer_m[0][2].diff(x)
        + u_outer_m[1][2].diff(y)
        + u_outer_m[2][2].diff(z)
    )
];
div_rho_inverse_m_outer_m = div_u_outer_m; # Exact for analytic expressions
