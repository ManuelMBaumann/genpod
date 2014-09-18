function eval = fun_ode45(t,y)

global Bred Cred nu Upod

eval = -0.5*Bred*(Upod*y).^2 - nu*Cred*y;

end