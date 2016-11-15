
close all;
h = -100:.001:100;
C = zeros(size(h));
K = zeros(size(h));
theta =zeros(size(h));

for i= 1:length(h)
    [C(i),K(i),theta(i)] = vanGenuchten(h(i),phi);   
end

plot(C, h, '.' )

