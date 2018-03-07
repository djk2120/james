clear
close all


thedir = '../data/mar6/';
x = dir(thedir);
files = cell(length(x)-2,1);
for i = 3:length(x)
    files(i-2) = {[thedir,x(i).name]};
end

files = files([5,6,1]);

offset = 10;
nx = length(files);
n  = length(ncread(files{1},'FCTR'))-offset;

a  = 1;
ct = 0;
month = zeros(1,48*365);
for x = eomday(2001,1:12)
ct = ct+1;
month(a:a+x*48-1) = ct;
a = a +x*48;
end
month = repmat(month,1,3);
month = month(1:n);
year  = 1:length(month);
year(year<48*365+1) = 1;
year(year>48*365*2) = 3;
year(year>3)        = 2;
tt    = repmat(1:48,1,365*3);
tt    = tt(1:n);

zs=[     0  , 0.0200  , 0.0600  , 0.1200  , 0.2000  , 0.3200  , 0.4800,...
    0.6800  , 0.9200  , 1.2000  , 1.5200  , 1.8800  , 2.2800  , 2.7200,...
    3.2600  , 3.9000  , 4.6400  , 5.4800  , 6.4200  , 7.4600  , 8.6000];

zr=[1.40e-2,2.73e-2,3.96e-2,5.02e-2,7.02e-2,...
    8.49e-2,9.36e-2,9.62e-2,9.36e-2,8.67e-2,...
    7.68e-2,6.54e-2,5.36e-2,4.67e-2,3.67e-2,...
    2.62e-2,1.71e-2,1.03e-2,5.70e-3,2.92e-3];

if ~exist('qsink','var')
h2osoi = zeros(20*nx,n);
qsink  = h2osoi;
ksr    = h2osoi;
et     = zeros(nx,n);
vwp    = zeros(4*nx,n);
for ff=1:nx
    tmp = ncread(files{ff},'H2OSOI');
    h2osoi((1:20)+(ff-1)*20,:) = tmp(1,:,1+offset:end);
    tmp = ncread(files{ff},'FCTR');
    et(ff,:) = tmp(1+offset:end);
    tmp = ncread(files{ff},'VEGWP');
    vwp((1:4)+(ff-1)*4,:) = tmp(1,1:4,1+offset:end);
    tmp = ncread(files{ff},'QROOTSINK');
    qsink((1:20)+(ff-1)*20,:) = tmp(1,:,1+offset:end);
    tmp = ncread(files{ff},'KSR');
    ksr((1:20)+(ff-1)*20,:) = tmp(1,:,1+offset:end);
end
end
%-----------------------------

rr = [1,0,0,0,0];

if rr(1)>0
figure
for i=[1,6]
    x = splitapply(@mean,et(i,:),month+(year-1)*12);
    if i==1
    plot(x,'k','LineWidth',1.5)
    else
        plot(x,'k:','LineWidth',1.5)
    end
    hold on
    ylim([0 150])
    xlim([0 36])
    set(gca,'ytick',0:50:150)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
    grid on
end

end

if rr(2)>0
figure
for i=1:nx
    zv = zs(2:end)-zs(1:end-1);
    x=zv(1:13)*h2osoi((1:13)+(i-1)*20,:)+h2osoi(14+(i-1)*20,:)*(3-zs(14));
    t= 2001+(1:length(x))/(48*365);
    plot(t,1000*x,'LineWidth',1.5)
    set(gca,'xtick',2001:2004)
    hold on
end
ylim([200,1000])
end

if rr(3)>0
    
    for i=1:8
    subplot(2,4,i)
    ix = year==2&month==11;
    t  = 0.5:0.5:24;
    x  = splitapply(@mean,vwp(1+(i-1)*4,ix),tt(ix));
    plot(t,x/101972)
    hold on
    x  = splitapply(@mean,vwp(4+(i-1)*4,ix),tt(ix));
    plot(t,x/101972)
    hold on
    set(gca,'xtick',0:6:24)
    xlim([0 24])
    ylim([-3.3,0])
    end
              
end

if rr(4)>0
    hr = qsink;
    for i=1:length(qsink(:,1))
        hr(i,:) = -qsink(i,:).*(qsink(i,:)<0);
    end
    
    hr2 = zeros(nx,n);
    for i=1:nx
        hr2(i,:) = sum(hr((1:20)+(i-1)*20,:));
    end
    
    for i=1:nx
        plot(1800*cumsum(hr2(i,:)))
        hold on
    end
    legend(num2str((1:nx)'))
    
end

if rr(5)>0
    for xx=1:2
    for ss=2:10
        targ = ksr(ss+(xx-1)*20,:);
        xmin = min(targ);
        xmax = max(targ);
        dx   = (xmax-xmin)/501;
        xv   = xmin+dx:dx:xmax-dx;
        out  = zeros(500,1);
        last = 0;
        for i=1:500
            a  = sum(targ<xv(i))/n;
            out(i) = a;
            last = a;
        end
        subplot(3,3,ss-1)
        plot(xv,out)
        hold on
    end
    end
end
    
    
    
    
    



    