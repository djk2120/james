function [ a ] = getvars( file , offset, ns, nosnow)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
temp  = double(ncread(file,'mcdate'));
n     =length(temp)-offset;
mcdate    = zeros(1,n);
mcdate(:) = temp(1:end-offset);

mdcur     = zeros(1,n);
mcsec     = zeros(1,n);
prec      = zeros(1,n);
rh        = zeros(1,n);
fsds      = zeros(1,n);
flds      = zeros(1,n);
pbot      = zeros(1,n);
zbot      = zeros(1,n);
tbot      = zeros(1,n);
windv     = zeros(1,n);

temp      = double(ncread(file,'mdcur'));
mdcur(:)  = 1+temp(1:end-offset);
temp      = double(ncread(file,'mcsec'));
mcsec(:)  = 1+temp(1:end-offset);
if nosnow
    temp = ncread(file,'RAIN');
else
    temp=ncread(file,'RAIN')+ncread(file,'SNOW');
end
prec(:)=temp(1,1+offset:end);
temp=ncread(file,'RH');
rh(:)=temp(1,1+offset:end);
temp=ncread(file,'FSDS');
fsds(:)=temp(1,1+offset:end);
temp=ncread(file,'FLDS');
flds(:)=temp(1,1+offset:end);
temp=ncread(file,'PBOT');
pbot(:)=temp(1,1+offset:end);
temp=ncread(file,'ZBOT');
zbot(:)=temp(1,1+offset:end);
temp=ncread(file,'TBOT');
tbot(:)=temp(1,1+offset:end);
temp=ncread(file,'WIND');
windv(:)=temp(1,1+offset:end);

tc = tbot - 273.15;
vpd = (1-rh/100).*(0.61078*exp((17.27*tc)./(tc+237.3)));


year  = floor(mcdate/10000);
month = floor((mcdate-10000*year)/100);
day   = mcdate-10000*year-100*month;




y1=1+int8(year-min(year));
y2=year;
y2(1:end-1)=y2(2:end);
idx=year~=y2;
endlist=[0,mdcur(idx)];
doy=mdcur-endlist(y1);


ysh=zeros(1,mdcur(end));
msh=zeros(1,mdcur(end));
daysh=zeros(1,mdcur(end));
doysh=zeros(1,mdcur(end));
for dd=1:max(mdcur)
ysh(dd)   = min(year(mdcur==dd));
msh(dd)   = min(month(mdcur==dd));
daysh(dd) = min(day(mdcur==dd));
doysh(dd) = min(doy(mdcur==dd));
end

ylist = unique(year);
dlist = unique(mdcur);
diurn = unique(mcsec);
nt    = length(diurn);

dailyp = zeros(1,mdcur(end));
todayp = 0*prec;
p30    = 0*prec-1;
p30sh  = zeros(1,max(mdcur))-1;
wind   = 30;
for dd=1:max(mdcur)-1
    idx         = mdcur==dd;
    dailyp(dd) = 1800*sum(prec((dd-1)*nt+(1:48)));
    todayp(idx) = dailyp(dd);
    if dd>wind
        p30(idx) = sum(dailyp(dd-wind:dd-1));
        p30sh(dd)= sum(dailyp(dd-wind:dd-1));
    end
end




assignin('caller','mcdate',mcdate);
assignin('caller','n'     ,n);
assignin('caller','mdcur' ,mdcur);
assignin('caller','mcsec' ,mcsec);
assignin('caller','prec'  ,prec);
assignin('caller','rh'    ,rh);
assignin('caller','fsds'  ,fsds);

assignin('caller','flds'  ,flds);
assignin('caller','wind'  ,windv);
assignin('caller','pbot'  ,pbot);
assignin('caller','zbot'  ,zbot);
assignin('caller','tbot'  ,tbot);

assignin('caller','year'  ,year);
assignin('caller','month' ,month);
assignin('caller','day'   ,day);
assignin('caller','doy'   ,doy);
assignin('caller','ylist' ,ylist);
assignin('caller','dlist' ,dlist);
assignin('caller','diurn' ,diurn);
assignin('caller','dailyp',dailyp);
assignin('caller','todayp',todayp);
assignin('caller','p30'   ,p30);
assignin('caller','p30sh' ,p30sh);
assignin('caller','ysh'   ,ysh);
assignin('caller','msh'   ,msh);
assignin('caller','daysh' ,daysh);
assignin('caller','doysh' ,doysh);
assignin('caller','vpd',vpd);


a=ncinfo(file);


end

