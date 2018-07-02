clear
close all

xx = 0.6;
x = rand(1000,1);
tt = 0;
mt = 0;
for i=1:1000
    if x(i)<xx
        tt = tt+1;
    else
        mt = max(mt,tt);
        tt = 0;
    end
end



mt