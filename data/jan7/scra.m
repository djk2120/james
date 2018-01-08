clear
close all

file = 'BR-CAX_I1PTCLM50_r270.clm2.h1.2001-01-01-00000.nc';


t = ncread(file,'mdcur');
t   = t(1:end-1)';



day = double(t)+1;
year = 0*day+1;
year(day>365) = year(day>365)+1;
day(day>365) = day(day>365)-365;
year(day>365) = year(day>365)+1;
day(day>365) = day(day>365)-365;

month=day;
mvals = [0,cumsum(eomday(2001,1:12))];
for i=1:12
    month(day>mvals(i)&day<=mvals(i+1))=i;
end
tt = repmat(1:48,1,1095);


tmp = ncread(file,'FCTR');
targ = 0*month;
targ(:) = tmp(2:end);

ix = year==3&month==10;
sum(ix)

out = splitapply(@mean,targ(ix),tt(ix));

plot(out)
