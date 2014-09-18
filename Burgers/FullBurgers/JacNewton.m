function J=JacNewton(y)

global M C nu dt

% J = (1/dt)*M + Ny(y) + nu*C;
J = (1/dt)*M + 0.5*( Ny(y) + nu*C );

end