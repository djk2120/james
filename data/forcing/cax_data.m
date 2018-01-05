clear
close all

file = 'clmforc.1x1pt_Br-CAX.nc.orig';
newf = 'clmforc.1x1pt_Br-CAX.nc';
offset = 0;

aa = ncinfo(file);



t = ncread(file,'time');
t   = t';

day = floor(t)+1;
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


%check the variables

varlist = {'TBOT','RH','WIND','FLDS','FSDS','PRECTmms'};

if 1==2
    for i=1:6
        
        tmp = ncread(file,varlist{i});
        targ = zeros(1,length(tmp));
        targ(:) = tmp(1,1,:);
        
        ix  = month==10&year==3;
        out = splitapply(@mean,targ(ix),tt(ix));
        
        subplot(2,3,i)
        plot(out)
        xlim([0 49])
        title(varlist{i})
    end
end


if 1==2
    i=5;
    tmp = ncread(file,varlist{i});
    targ = zeros(1,length(tmp));
    targ(:) = tmp(1,1,:);
    
    d1 = 274;
    
    for i=1:4
        subplot(3,4,i)
        plot(targ(year==1&day==(d1+10*i)))
        subplot(3,4,i+4)
        plot(targ(year==2&day==(d1+10*i)))
        subplot(3,4,i+8)
        plot(targ(year==3&day==(d1+10*i)))
    end
end


%starts day 979, persists through end of record

ix = [];
a = 1:2:48;
b = 2:2:48;
for i=1:24
    ix = [ix;b(i);a(i)];
end


%starts being bad
%tbot = 954
%rh   = 954
%wind = the whole time?
%fsds = 979


if 1==2
    i=3;
    tmp = ncread(file,varlist{i});
    targ = zeros(1,length(tmp));
    targ(:) = tmp(1,1,:);
    targ2 = targ;
    
    %for dd=979:1095
    xdk = figure;
    xdk.Units = 'inches';
    xdk.Position = [.2639  , 2.1528 , 16.1389  ,  7.6250];
    ct = 0;
    for dd=900:970
        ii = (1:48)+(dd-1)*48;
        x  = targ(ii);
        %fsds2(ii) = x(ix);
        
        ct = ct+1;
        if ct==7
            ct=1;
        end
        
        subplot(2,6,ct)
        plot(x)
        title(num2str(dd))
        subplot(2,6,ct+6)
        plot(x(ix))
        pause(0.5)
        
    end
end



if 1==1
    xdk = figure;
    xdk.Units = 'inches';
    xdk.Position = [.2639  , 2.1528 , 16.1389  ,  7.6250];
    
    i=3;
    tmp = ncread(file,varlist{i});
    targ = zeros(1,length(tmp));
    targ(:) = tmp(1,1,:);
    
    for dd=1:1095
        ii = (1:48)+(dd-1)*48;
        x  = targ(ii);
        targ2(ii) = x(ix);
    end
      
    ct = 0;
    for yy=1
        for mm=1:12
            ct = ct+1;
            if ct==7
                ct=1;
            end
            ix  = year==yy&month==mm;
            out = splitapply(@mean,targ(ix),tt(ix));
            subplot(2,6,ct)
            plot(out)
            title([num2str(mm),'-',num2str(yy+2000)])
            out = splitapply(@mean,targ2(ix),tt(ix));
            subplot(2,6,ct+6)
            plot(out)
            pause(0.5)
        end
    end
    
    
    
end
