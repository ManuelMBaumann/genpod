function Yapprox = get_red_model(Yfull, y0, flag)

global M B C Nt N k1 k2 ns dt T plotOn Upod Cred Bred Fred


if flag == 1
    addpath('./POD');
    %disp('POD')
elseif flag == 2
    addpath('./POD');
    addpath('./genPOD');
    %disp('genPOD')
else
    addpath('./PODDEIM');
    %disp('POD-DEIM')
end


%% Setup POD
if flag == 1
    
    %R = chol(M);
    [Uy,Sy,~]=svd(Yfull);
    %Upod=(R')\Uy(:,1:k1);
    Upod = Uy(:,1:k1);
    
    if plotOn ==1
       figure(55)
       index=1:min(N-1,Nt);
       semilogy(index,diag(Sy(index,index)),'x'); %investigate basis dimensions
    end


    % Pre-compute matrices
    Cred  = Upod'*C*Upod;
    Bred  = Upod'*B;
    Mred  = Upod'*M*Upod;
    
    invM = inv(Mred);
    Cred = invM*Cred;
    Bred = invM*Bred;
    
    % obtain reduced variables
    y0_red = Upod'*y0;
    
elseif flag == 2
    
    % Compute mass matrix Ypl = < y_i , \nu_j > and My = < \nu_i , \nu_j >
    [Ypl,My,~] = massPL(N-1,ns,0,T,Yfull);
    
    %Compute eigs of R
    %[Z,D]=eig(inv(My));
    %M12=Z*sqrt(D)*Z';
    R = chol(My);
    
    %norm(inv(R)-M12)
    
    [Upod,~,~]=svd(Ypl*inv(R));
    
    Upod = Upod(:,1:k1);
      
    % Pre-compute matrices
    Cred  = Upod'*C*Upod;
    Bred  = Upod'*B;
    
    Mred = Upod'*M*Upod;
    
    invM = inv(Mred);
    Cred = invM*Cred;
    Bred = invM*Bred;
    
    y0_red = Upod'*y0;
    
else

    %% Setup POD-DEIM

    % SVD of snapshot matrices + DEIM
    R = chol(M);
    [Uy,Sy,~]=svd(R*Yfull);
    Upod=R\Uy(:,1:k1);
    
    [Uf,Sf,~]=svd(Yfull.^2);
    Udeim=Uf(:,1:k2);
    [~,ind]=deim(Udeim);

    if plotOn == 1
       figure(55)
       index=1:min(N-1,Nt);
       semilogy(index,diag(Sy(index,index)),'x',index,diag(Sf(index,index)),'rx'); %investigate basis dimensions
    end
    
    % Pre-compute matrices
    Bred  = Upod'*B*Udeim/Udeim(ind,:);
    Cred  = Upod'*C*Upod;
    Fred  = Upod(ind,:);
    
    % obtain reduced variables
    y0_red = Upod'*M*y0; 
end
    
%% Solve red system
% tStart = tic;
% Y_red = Burgers_Reduced(y0_red);
% tElapsed_red = toc(tStart);  

[~, Y_red] = ode45(@fun_ode45,0:dt:T,y0_red);    
Yapprox = Upod*Y_red';

end






