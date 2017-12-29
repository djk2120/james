function [ good,le,ta,goodt,vpd,goodv,gpp,goodg ] = get_obs(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%file   = '/Users/kennedy/work/gmd/data/obs/AMF_USMe2_2004_L2_WG_V005.nc';
%a      = ncinfo(file);

f1     = '../data/conr207calwspinsp_US-Me2_I1PTCLM50.clm2.h1.2002-01-01-00000.nc';
offset = 17;

if ~exist('mcdate','var')
    % get the vars that don't vary by file
    a = getvars( f1 , offset, 11);
end

%if ~exist('fctr','var')
    % get the vars that DO vary by file
   % varlist={'FCTR','TA'};
    %           1      2      3      4       5     6     7
  %  vard=ones(length(varlist),1);
 %   x = getmore( {f1},offset,n,varlist,vard );  
%end


le = zeros(1,160000);
ta = zeros(1,160000);
vpd= zeros(1,160000);
ct=1;
for yy = 2002:2010
    file   = ['/Users/kennedy/work/gmd/data/obs/AMF_USMe2_',num2str(yy),'_L2_WG_V005.nc'];
    
    temp                      = ncread(file,'LE');
    le(ct:ct+length(temp)-1)  = temp(:);
    temp                      = ncread(file,'TA');
    ta(ct:ct+length(temp)-1)  = temp(:);
    temp                      = ncread(file,'VPD');
    vpd(ct:ct+length(temp)-1) = temp(:);
    ct=ct+length(temp);
    
end


gpp= zeros(1,160000)-333;
ct=1;
for yy = 2002:2007
    file   = ['/Users/kennedy/work/gmd/data/obs/AMF_USMe2_',num2str(yy),'_SYN_WG_V000.nc'];
    temp                      = ncread(file,'GPP');
    gpp(ct:ct+length(temp)-1)  = temp(:);
    ct=ct+length(temp);
    
end



le      = le(1:length(mcdate));
good    = le>-9000;

ta      = ta(1:length(mcdate));
goodt   = ta>-100;

vpd     = vpd(1:length(mcdate));
goodv   = vpd>-100;

gpp     = gpp(1:length(mcdate));
goodg   = gpp>-100;

end

