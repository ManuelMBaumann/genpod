function [errM, gerrM] = MOR_driver(ns_input,k1_input)

global nu tol_newton max_newton h dt M B C Nt N k1 k2 ns T ...
    time space plotOn 

ns = ns_input;
k1 = k1_input;

addpath('./FullBurgers');

%% Params
N  = 100;                % # grid point in space
Nt = 500;               % # grid point in time
L  = 1;                 % x \in [0,L]
T  = 1;                 % t \in [0,T]

x = linspace(0,L,N+1);           % Spacial discretization
t = linspace(0,T,Nt);            % time discretization
[time, space] = meshgrid(t',x);  % for 3D plotting
h = L/N;                         % constant step size
dt= T/(Nt-1);                    % constant time step

nu    = 0.01;           % viscosity parameter

tol_newton = 1e-12;     % for Newton's method
%tol_newton_red = 1e-12;
max_newton = 40;        % max. number of inner Newton iterations
%max_newton_red = 40;

%ns = 17;
%k1 = 17;                    % # POD basis of y
% k2 = 25;                  % # DEIM basis

ss = Nt/ns;
ss_vec = floor(1:ss:Nt);
plotOn = 0;             % plotting and output on (1) or off(0)


%% Initialize the state y_0
y0 = [ones(floor((N-1)/2),1);zeros(ceil((N-1)/2),1)];
%y0 = sin(linspace(0,L,N-1)*2*pi/L)';

% Pre-compute stiffness and mass matrix for full-order Burgers eqn
M = (h/6)*gallery('tridiag',ones(N-2,1),4*ones(N-1,1),ones(N-2,1));
B = gallery('tridiag',-0.5*ones(N-2,1),zeros(N-1,1),0.5*ones(N-2,1));
C = (1/h)*gallery('tridiag',-ones(N-2,1),2*ones(N-1,1),-ones(N-2,1));


%% Solve Full-order system
Yfull = Burgers(y0);
%Yfull = Yfull + 0.001*rand(N-1,Nt);

Yapprox1 = get_red_model(Yfull(:,ss_vec), y0, 1);
Yapprox2 = get_red_model(Yfull, y0, 2);

%% Plot full-order system and reduced model in one figure
figure(1)
subplot(1,3,1);
mesh(time, space, [zeros(1,Nt);Yfull;zeros(1,Nt)]); %add D-BC to Yfull
xlabel('t')
ylabel('x')
zlabel('y(t,x)'); %zlim([0 1.1])
title(['Full Order System, n = ' num2str(N)])
subplot(1,3,2);
mesh(time, space, [zeros(1,Nt);Yapprox1;zeros(1,Nt)]); %add D-BC to Yapprox
xlabel('t')
ylabel('x')
zlabel('y_k(t,x)'); %zlim([0 1.1])
title(['# POD basis = ' num2str(k1) '     # DEIM basis = ' num2str(k2)])
subplot(1,3,3);
mesh(time, space, [zeros(1,Nt);Yapprox2;zeros(1,Nt)]); %add D-BC to Yapprox
xlabel('t')
ylabel('x')
zlabel('y_k(t,x)'); %zlim([0 1.1])
title(['# genPOD basis = ' num2str(k1) '     # DEIM basis = ' num2str(k2)])


%% Error plot
figure(99)
mesh(space, time, [zeros(1,Nt);Yfull;zeros(1,Nt)]); %add D-BC to Yfull
xlabel('x')
ylabel('t')
title('Full Order Soultion', 'FontSize', 16);
shading interp
colorbar
axis equal
view(0,90)

figure(100)
%subplot(1,2,1);
errplotPOD = abs([zeros(1,Nt);Yapprox1;zeros(1,Nt)] - ...
    [zeros(1,Nt);Yfull;zeros(1,Nt)] );
surf(space, time, errplotPOD); %add D-BC to Yapprox
xlabel('x')
ylabel('t')
title(['POD error'])
shading interp
cmin = min(min(errplotPOD)); cmax = max(max(errplotPOD));
caxis([cmin cmax/5])
colorbar
colormap cool
axis equal
view(0,90)

figure(101)
%subplot(1,2,2);
errplotPODgen = abs([zeros(1,Nt);Yapprox2;zeros(1,Nt)] - ...
    [zeros(1,Nt);Yfull;zeros(1,Nt)]);
surf(space, time, errplotPODgen); %add D-BC to Yapprox
xlabel('x')
ylabel('t')
zlabel('error');
title(['genPOD error'])
shading interp
caxis([cmin cmax/5])
colormap cool
axis equal
colorbar
view(0,90)

%% errors
errM  = Mnorm( Yfull, Yapprox1, M );
gerrM = Mnorm( Yfull, Yapprox2, M );

disp(['genpod vs. pod error: ', num2str(gerrM), ' ',  num2str(errM)])
disp(['genpod vs. pod relative change: ', num2str(100*(gerrM-errM)/gerrM)])


%% remove paths for next run
warning off
rmpath('./FullBurgers');
rmpath('./POD');
rmpath('./genPOD');
rmpath('./PODDEIM');

end