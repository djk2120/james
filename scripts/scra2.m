clear
close all

f='../data/k4g7.nc';

pftname = ncread(f,'pftname');
pftname(:,5)'

s = ncread(f,'smpso');
s(5)/101972

s = ncread(f,'smpsc');
s(5)/101972