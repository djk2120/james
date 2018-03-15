clear
close all

p = -0.01:-0.01:-5;

p50 = -2.5;

f1 = @(x,a) 1./(1+(x/p50).^a);
f2 = @(x,a) 2.^-((x/p50).^a);

plot(p,f1(p,6))
hold on
plot(p,f2(p,5))