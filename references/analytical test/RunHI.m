%Call HortonImpervious
tss = [60,90,120,150,300,600];
for i = 1:6

So = 0.05;              % rise/run
n = 0.03;               % dimensionless
P = 30/3600/100;        %m/s
L = 100;                % m
%ts = 0.5*3600;          % s
ts=tss(i);
nfall = 1000;           % Number of points solved in falling limb

[Q,t]=HortonImpervious(So,n,P,L,ts,nfall);

plot(t,Q); hold on
end