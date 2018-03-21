clear
close all


p = -0.1:-0.1:-5;

p50 = -2.5;
ck  = 4.95;
a   = 6


f1 = 2.^-((p/p50).^ck);
plot(p,f1)

hold on
f2 = 1./(1+(p/p50).^a)
plot(p,f2)