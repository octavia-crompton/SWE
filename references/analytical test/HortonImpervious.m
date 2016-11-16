function [Q,t]=HortonImpervious(So,n,P,L,ts,nfall)
%%%%%%%%%%%%%%%%
% This code solves for a Hortonian Hillslope Hydrograph
% Under the assumption of zero infiltration
% The inputs are the slope angle -- So
% Slope roughness (Manning n) -- n
% Precipitation rate (m/s) -- P
% Slope length (m) -- L
% Storm duration (s) -- ts
% Because the falling limb is solved implicitly, we can also specify the
% Number of solution points.  Increasing this from its default may better
% Resolve the falling limb.

% It is assumed that the time of ponding is t=0
% and that the storm begins at t=0

Kr = So^(1/2)*1/n;      % Simplifies constants to a single kinematic roughness

te = (L./(Kr * P.^(2/3))).^(3/5);

if ts>te
    
    trise = linspace(0,te,100);     % We'll get 100 flow points on the rising limb
    Qrise = Kr.*(P.^(5/3)).*(trise.^(5/3)); % And generate flow for each of them
    
    teq = linspace(te,ts,5);        % We'll generate 5 points under equilibrium conditions
    Qeq = P.*L.*ones(size(teq));    % And generate flow for each of those
    
    Qfall = linspace(Qeq(end),0,nfall);    % We'll get nfall flow points on the falling limb
    tfall = ts+(L-Qfall./P)./(5/3*Kr.^(3/2).*Qfall.^(2/5));
    
else
    th = (L-Kr./P*(P.*ts).^(5/3))/(5/3 * Kr * (P.*ts).^(2/3));
    
    trise=linspace(0,ts,100);
    Qrise = Kr.*(P.^(5/3)).*(trise.^(5/3));
    
    teq = linspace(ts,ts+th,5);
    Qeq = Kr.*(P.*ts).^(5/3)*ones(size(teq));
    
    Qfall = linspace(Qeq(5), 0,nfall);
    tfall = ts+th+(L-Qfall./P)./(5/3*Kr.^(3/2).*Qfall.^(2/5));
end

Q=[Qrise,Qeq,Qfall];
t=[trise,teq,tfall];

