function Y_red = Burgers_Reduced(y0_red)

global tol_newton_red Nt k1 max_newton_red

rhs_red = zeros(k1,1); % zero rhs in Burgers eqn

Y_red      = zeros(k1,Nt);
Y_red(:,1) = y0_red;

for tt=2:Nt %time loop
    y_tmp1=Y_red(:,tt-1); %nice guess  
    err=1; 
    iter = 1;
    while (err > tol_newton_red && iter < max_newton_red)%inner Newton method
       f_red=fNewton_POD(Y_red(:,tt-1),y_tmp1,rhs_red); 
       J_red=JacNewton_POD(y_tmp1);

       ytilde_red=J_red\-f_red;
       y_tmp2=y_tmp1+ytilde_red;

       err=norm(y_tmp2-y_tmp1)/norm(y_tmp2);
       y_tmp1=y_tmp2;  
       
       iter = iter+1;
    end
    Y_red(:,tt)=y_tmp2;
end


end

