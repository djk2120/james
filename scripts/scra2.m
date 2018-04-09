clear
close all

n = 999999;
a=randn(n,1);

kk = 50000+a*10000;

p10 = 37184;

revs = zeros(n,1);
ix   = kk>p10;
revs(ix)  = p10*0.2+(kk(ix)-p10)*0.05;
revs(~ix) = kk(~ix)*0.2;

mean(revs)