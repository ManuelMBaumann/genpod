function eval = fNewton(Yold,y,f)

global h M B C nu dt

% eval = (1/dt)*M*y - (1/dt)*M*Yold + 0.5*B*y.^2 + nu*C*y - h*f;
eval = (1/dt)*M*y - (1/dt)*M*Yold - 0.5 * ( -0.5*B*y.^2 - nu*C*y + h*f ...
    - 0.5*B*Yold.^2 - nu*C*Yold + h*f);

end