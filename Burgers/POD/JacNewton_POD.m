function J_red=JacNewton_POD(y)

global Cred nu dt k1

J_red = (1/dt)*eye(k1) + Ntildey(y) + nu*Cred;

end