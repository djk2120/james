clear
close all

file = 'clmforc.1x1pt_Br-CAX.nc';
offset = 0;

a = ncinfo(file);


t = ncread(file,'time');
t   = t';

tmp = ncread(file,'FSDS');
fsds = zeros(1,length(tmp));
fsds(:) = tmp(1,1,:);


%starts day 979, persists through end of record
fsds2 = fsds;
ix = [];
a = 1:2:48;
b = 2:2:48;
for i=1:24
ix = [ix;b(i);a(i)];
end

for dd=979:1095
    ii = (1:48)+(dd-1)*48;
    x  = fsds(ii);
    fsds2(ii) = x(ix);
end

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
out = zeros(48,1);
for yy=1:3
    for mm=11
        ix = year==yy&month==mm;
        out = splitapply(@mean,fsds2(ix),tt(ix));

        
        hold on
        plot(out)
        title([num2str(mm),'-',num2str(2000+yy)])
        
    end
end
    
figure
a = sum(eomday(2001,1:10));
for dd = 1:30
    ix = year==3&month==11&day==(a+dd);
    subplot(1,2,1)
    plot(fsds(ix))
    ylim([0 1000])
    subplot(1,2,2)
    plot(fsds2(ix))
    title(dd)
    ylim([0 1000])
    pause(0.8)
end
    





