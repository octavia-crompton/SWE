%  Note: Axes are switched  (correct orientation)
%  No surface flux 
% Configured to set initial conditions to the output of the most recent run

% DESCRIPTION: This is a 1-D Richards equation solver coded up for GEOS   %
% 697: Fundamentals of Simulation Modeling in the Hydrologic Sciences at  %
% Boise State University for the Fall 2010 semester. It is intended for   % 
% student use only.                                                       %
%                                                                         %
% The formulation is based on a backward Euler implementation of the      %
% mixed theta-head formulation of the 1-D Richards Equation, with a       %
% modified Picard iteration finite differencing solution scheme. This     %
% implementation follows exactly the algorithm outlined by Celia et al.   %
% [1990] in Water Resources Research.                                     %
%                                                                         %
% The soil water retention parameterization used is that of van Genuchten %
% [1980], and this code requires the corresponding MATLAB function        %
% vanGenuchten.m.                                                         %
%                                                                         %
% References:                                                             %
%                                                                         %
% Celia, M.A., Bouloutas, E.T. and Zarba, R.L., 1990. A general mass-     %
%     conservative numerical solution for the unsaturated flow equation.  %
%     Water Resources Research 26(7), 1483?1496.                          %
%                                                                         %
% van Genuchten, M.Th., 1980. A closed-form equation for predicting the   %
%     hydraulic conductivity of unsaturated soils. Soil Science Society   %
%     of America Journal 44, 892?898.                                     %
%                                                                         %
%=========================================================================%% Richard's Solver
% clear all;
% close all;
% Stopping tolerance [cm]
stop_tol = 0.01;
% Define van Genuchten parameters
alpha   = 0.0335;
theta_S = 0.368;
theta_R = 0.102;
lambda  = 0.25;
n       = lambda + 1;
m       = lambda/n;
Ksat    = 0.09;

phi(1) = alpha;
phi(2) = theta_S;
phi(3) = theta_R;
phi(4) = n;
phi(5) = m;
phi(6) = Ksat;

% Make mesh variables
dz = 1; % [cm]
zmin = 0; % [cm]
zmax = 100; % [cm]
z    = (zmin:dz:zmax);
%% 
nz   =  length(z);

% Define time variables
dt = 1800; % [s]
tmin = 0;  % [s]
tmax = 28000; % [s]
t = (tmin:dt:tmax);
nt = length(t); 

% Define boundary conditions at top and bottom
htop = -100;
hbottom = -1000;
hinit = hnp1m;
% hinit(1) = htop;
% hinit(nz) = hbottom;

% Define matrices that we'll need in solution
BottomBoundCon = 1; % 0 = constant head, 1 = free drainage
DeltaPlus  = diag(-ones(nz,1)) + diag(ones(nz-1,1),1);
DeltaPlus(1,:) = 0;
DeltaPlus(nz,:) = 0;

DeltaMinus = diag(ones(nz,1)) + diag(-ones(nz-1,1),-1);
DeltaMinus(1,:) = 0;
DeltaMinus(nz,:) = 0;

MPlus  = diag(ones(nz,1)) + diag(ones(nz-1,1),1);
MPlus(1,1) = 2;
MPlus(1,2:nz) = 0;
MPlus(nz,nz) = 2;
MPlus(nz,1:(nz-1)) = 0;

MMinus = diag(ones(nz,1)) + diag(ones(nz-1,1),-1);
MMinus(1,1) = 2;
MMinus(1,2:nz) = 0;
MMinus(nz,nz) = 2;
MMinus(nz,1:(nz-1)) = 0; % Define initial conditions

% Define a storage container to store the pressure heads and soil moistures
[Ctop,Ktop,thetatop] = vanGenuchten(htop,phi);
[Cbottom,Kbottom,thetabottom] = vanGenuchten(hbottom,phi);
[Cinit,Kinit,thetainit] = vanGenuchten(hinit,phi);
H = zeros(nz,nt);
K = zeros(nz,nt);
H(:,1) = hinit;
K(:, 1) = Kinit;
THETA = zeros(nz,nt);
THETA(:,1) = thetainit;

figure(1); hold on;
plot(hinit, -flipud(z),'r');
xlabel('Pressure head [cm]');
ylabel('Depth [cm]');
title('Original pressure head');

figure(2); hold on;
plot(thetainit,-flipud(z),'r');
xlabel('Soil moisture [cm^3/cm^3]');
ylabel('Depth [cm]');
title('Original soil moisture');

% Define the container for an iteration counter
iterations = zeros(nt-1,1);
tic
for i=2:nt % Initialize the Picard iteration solver
hn  = H(:,i-1);
hnp1m = hn; 
thetan = THETA(:,i-1); % Enter the iteration step

% Define a dummy stopping variable
stop_flag = 0; 
% Define an iteration counter
niter = 0; 
while(stop_flag==0) % Get C,K,theta
    [cnp1m,knp1m,thetanp1m] = vanGenuchten(hnp1m,phi); 
    % 1. Compute the individual elements of the A matrix for LHS
    C = diag(cnp1m); 
    kbarplus = (1/2)*MPlus*knp1m;
    Kbarplus = diag(kbarplus); 
    kbarminus = (1/2)*MMinus*knp1m;
    Kbarminus = diag(kbarminus);
    A = (1/dt)*C - 1/((dz)^2)*(Kbarplus*DeltaPlus ...
        - Kbarminus*DeltaMinus); 
    % 2. Compute the residual of MPFD (RHS)
    R_MPFD = (1/(dz^2))*(Kbarplus*DeltaPlus*hnp1m - Kbarminus*DeltaMinus*hnp1m) ...
        + (1/dz)*(kbarplus - kbarminus) - (1/dt)*(thetanp1m - thetan); 
    % 3. Compute deltam for iteration level m+1
    deltam = pinv(A)*R_MPFD; % Increment iteration counter and display number of iterations
    niter = niter + 1;
    %disp(['Completed niter = ',int2str(niter), ' iterations.']); % Compute error metric
    if(max(abs(deltam(2:(nz-1))))<stop_tol)
        stop_flag = 1;
        hnp1mp1 = hnp1m + deltam; % Force boundary conditions  
        hnp1mp1(1) = hnp1mp1(2)-dz;
        if(BottomBoundCon==0)
            hnp1mp1(nz) = hbottom;
        elseif(BottomBoundCon==1)
            hnp1mp1(nz) = hnp1m(nz-1);
        end
        [cnp1m,knp1m,thetanp1m] = vanGenuchten(hnp1mp1,phi);
        thetanp1mp1 = thetanp1m;
    else
        hnp1mp1 = hnp1m + deltam;
        hnp1m = hnp1mp1; % Force boundary conditions
        hnp1m(1) = hnp1m(2) - dz;
        if (BottomBoundCon==0)
            hnp1m(nz) = hbottom;
        elseif(BottomBoundCon==1)
            hnp1m(nz) = hnp1m(nz-1);
        end
    end 
end

THETA(:,i) = thetanp1mp1; 
H(:,i) = hnp1mp1;
K(:,i)= knp1m;
if(mod(i,4)==0)
    figure(1); hold on;
    plot(hnp1mp1,-flipud(z));
    xlabel('Pressure head [cm]');
    ylabel('Depth [cm]');
    lgd = legend(num2str(0), num2str(t(4)),  num2str(t(8)), ...
          num2str(t(12)),num2str(t(16)), 'Location','northwest');
    figure(2); hold on;
    plot(thetanp1mp1,-flipud(z));
    xlabel('Soil moisture [cm^3/cm^3]');
    ylabel('Depth [cm]');
    %pause(1)
end 
 % Save number of iterations
iterations(i-1) = niter;
end
toc

Kp5 = (K(1, :) + K(2, :))/2.;
Kn5 = (K(end-1, :) + K(end-2, :))/2.; 
fluxin =  -  Kp5.*((H(2, :) - H(1, :))/dz + 1.)*dt ;
fluxout =   Kn5.*((H(end,:) - H(end-1, :))/dz + 1.)*dt ;
newmass = sum(THETA(:, 2:end) - THETA(:, 1:end-1) , 1);


% plot(fluxin(2:end) + fluxout(2:end) + newmass)
% plot(hnp1mp1,-flipud(z),'r');
% xlabel('Pressure head [cm]');
% ylabel('Depth [cm]');
% 
% figure(2); hold on;
% plot(thetanp1mp1,-flipud(z),'r');
% xlabel('Soil moisture [cm^3/cm^3]');
% ylabel('Depth [cm]');
% plot(hnp1mp1,-flipud(z),'r');
% xlabel('Pressure head [cm]');
% ylabel('Depth [cm]');
% 
% figure(2); hold on;
% plot(thetanp1mp1,-flipud(z),'r');
% xlabel('Soil moisture [cm^3/cm^3]');
% ylabel('Depth [cm]');

% figure(3); hold on;
% imagesc(t,flipud(z),H);
% shading interp;
% xlabel('Time [s]');
% ylabel('Depth [cm]');
% xlim([tmin tmax]);
% ylim([zmin zmax]);
% cmaph = colormap;
% colormap(flipud(cmaph));
% colorbar;
% title('Pressure head [cm]');
% 
% figure(4); hold on;
% imagesc(t,flipud(z),THETA);
% shading interp;
% xlabel('Time [s]');
% ylabel('Depth [cm]');
% xlim([tmin tmax]);
% ylim([zmin zmax]);
% cmaph = colormap;
% colormap(flipud(cmaph));
% colorbar;
% title('Soil moisture [cm^3cm^{-3}]');
% 
% figure(5);
% plot(t(2:nt),iterations,'o-');
% xlabel('Time [s]');
% ylabel('Number of iterations');save('Richards_outputs.mat','H','THETA');
