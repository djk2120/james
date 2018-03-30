clear
close all


f = '../data/k4g7.nc';
f = '../data/clm5_params.c171117.nc';

aa= ncinfo(f);

pft = ncread(f,'pftname')';


kmax = ncread(f,'kmax');
kmax(5,:)
kmax(7,:)

krmax = ncread(f,'krmax');
krmax(5,:)
krmax(7,:)

p50 = ncread(f,'psi50');
p50(5,:)
p50(7,:)


g1 = ncread(f,'medlynslope');
g1(5,:)
g1(7,:)


g0 = ncread(f,'medlynintercept');
g0(5,:)
g0(7,:)