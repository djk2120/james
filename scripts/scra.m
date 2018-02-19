clear
close all

f1 = ['/Users/kennedy/Desktop/james/data/feb14/',...
'b5_k15e-3.nc'];
f2 = ['/Users/kennedy/Desktop/james/data/jan24/',...
    'BR-CAX_I1PTCLM50_r270.clm2.h1.2001-01-01-00000.nc'];
f3 = ['/Users/kennedy/Desktop/james/data/feb14/',...
'b7_k15e-3.nc'];


zs=[     0  , 0.0200  , 0.0600  , 0.1200  , 0.2000  , 0.3200  , 0.4800,...
    0.6800  , 0.9200  , 1.2000  , 1.5200  , 1.8800  , 2.2800  , 2.7200,...
    3.2600  , 3.9000  , 4.6400  , 5.4800  , 6.4200  , 7.4600  , 8.6000];

zv = zs(2:end)-zs(1:end-1);

tmp = ncread(f2,'H2OSOI');
h2osoi = zeros(60,length(tmp)-10);
h2osoi(1:20,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f1,'H2OSOI');
h2osoi(21:40,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f3,'H2OSOI');
h2osoi(41:60,:) = tmp(1,1:20,1+10:end);


tmp = ncread(f2,'SMP');
smp = zeros(60,length(tmp)-10);
smp(1:20,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f1,'SMP');
smp(21:40,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f3,'SMP');
smp(41:60,:) = tmp(1,1:20,1+10:end);

tmp = ncread(f2,'QROOTSINK');
qrootsink = zeros(60,length(tmp)-10);
qrootsink(1:20,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f1,'QROOTSINK');
qrootsink(21:40,:) = tmp(1,1:20,1+10:end);
tmp = ncread(f3,'QROOTSINK');
qrootsink(41:60,:) = tmp(1,1:20,1+10:end);



xdk = figure;
i = 1;
x=zv(1:13)*h2osoi((1:13)+(i-1)*20,:)+h2osoi(14+(i-1)*20,:)*(3-zs(14));
i = 2;
y=zv(1:13)*h2osoi((1:13)+(i-1)*20,:)+h2osoi(14+(i-1)*20,:)*(3-zs(14));
i = 3;
z=zv(1:13)*h2osoi((1:13)+(i-1)*20,:)+h2osoi(14+(i-1)*20,:)*(3-zs(14));

plot(x,'LineWidth',1.2)
hold on
plot(y,'LineWidth',1.2)
plot(z,'LineWidth',1.2)
set(gca,'xtick',0:48*365:3*48*365)
set(gca,'xticklabel',0:48*365:3*48*365)
title('Water content top 3m')
ylabel('mm h2o')
xlabel('year')
xlim([0,3*48*365])
set(gca,'xticklabel',2001:2004)

if 1==2
xdk.Units = 'inches';
xdk.Position = [2,2,4,3];
xdk.PaperSize = [4,3];
xdk.PaperPosition = [0,0,4,3];

print(xdk,'../figs/top3m','-dpdf')
end

if 1==1
figure
for i=1:480
barh(-1:-1:-20,h2osoi(41:60,i))
xlim([0,0.5])
pause(0.1)

end
end


tmp = ncread(f2,'FCTR')+ ncread(f2,'FCEV') + ncread(f2,'FGEV');
et = zeros(3,length(tmp)-10);
et(1,:) = tmp(1+10:end);
tmp = ncread(f1,'FCTR')+ ncread(f1,'FCEV') + ncread(f1,'FGEV');
et(2,:) = tmp(1+10:end);
tmp = ncread(f3,'FCTR')+ ncread(f1,'FCEV') + ncread(f1,'FGEV');
et(3,:) = tmp(1+10:end);

dd = repmat(1:365*3,48,1);
dd = dd(:);
dd = dd(1:52551);
etperday = 1800*4e-7*splitapply(@sum,et',dd);



hr = zeros(1,length(qrootsink));
ss = 40;
for i=1:20
    hr = hr+qrootsink(ss+i,:).*(qrootsink(ss+i,:)<0);
end
figure
plot(cumsum(hr))

    
 xf = 1000/101972;  %converts mm to kPa


t = 0.08:0.005:0.42;
s = t/0.42;
p = 500*s.^-4;

figure
%plot(t,-p)
if 1==1
hold on
for i=1:20
plot(h2osoi(40+i,24:48:end),smp(40+i,24:48:end)/101972,'.')
end
end





