clear
close all


t = 0.12:0.01:0.42;
tsat = 0.42;
psat = -100;
ksat = .015;
b    = 6;


p = psat*(t/tsat).^(-b);

%plot(t,log10(-p))


k = ksat*(t/tsat).^(2*b+3);
plot(t,log10(k))