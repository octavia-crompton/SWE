
close all;
h = -10:.1:10;
C = zeros(size(h));
K = zeros(size(h));
theta =zeros(size(h));

for i= 1:length(h)
    [C(i),K(i),theta(i)] = vanGenuchten(h(i),phi);   
end

plot(h, theta)
