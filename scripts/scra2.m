clear
close all

zs=[     0  , 0.0200  , 0.0600  , 0.1200  , 0.2000  , 0.3200  , 0.4800,...
    0.6800  , 0.9200  , 1.2000  , 1.5200  , 1.8800  , 2.2800  , 2.7200,...
    3.2600  , 3.9000  , 4.6400  , 5.4800  , 6.4200  , 7.4600  , 8.6000];
z = zs(2:end)-zs(1:end-1);


f = '../data/mar6/BR-CAX_I1PTCLM50_r270v3_phs_amb.clm2.h1.2001-01-01-00000.nc';

tmp = ncread(f,'H2OSOI');

s = zeros(20,365*3*48);
s(:,:) = tmp(1,:,1:365*3*48);

tmp = ncread(f,'FCTR');
t = tmp(1:365*3*48);

dd = repmat(1:365*3,48,1);
dd = [dd(:)]';

%top3m = [z(1:8),0.66]*s(1:9,:);
top = splitapply(@mean,[z(1:13),0.28]*s(1:14,:),dd);
bot = splitapply(@mean,[0.26,z(15:end)]*s(14:end,:),dd);

xdk = figure;

plot(top)
hold on
plot(bot)

ylim([0,max(2.5,max(bot))])
xlim([0,1095])

xdk.Units = 'inches';
xdk.Position = [2,2,7,4];
xdk.PaperSize = [7,4];
xdk.PaperPosition = [0,0,7,4];
print(xdk,'wat','-dpdf')

xdk = figure;
plot(splitapply(@mean,t,dd),'bx')
ylim([0,160])

xdk.Units = 'inches';
xdk.Position = [2,2,7,4];
xdk.PaperSize = [7,4];
xdk.PaperPosition = [0,0,7,4];
print(xdk,'t','-dpdf')

