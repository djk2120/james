close all

dir = '../data/mar6/';

files = {...
%    [dir,'BR-CAX_I1PTCLM50_r270.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r270_on_tfe.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r270_off_amb.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r270_off_tfe.clm2.h1.2001-01-01-00000.nc']};
[dir,'BR-CAX_I1PTCLM50_r270v3_phs_amb.clm2.h1.2001-01-01-00000.nc'];...
[dir,'BR-CAX_I1PTCLM50_r270v3_phs_tfe_60.clm2.h1.2001-01-01-00000.nc'];...
%[dir,'alt.nc'];...
[dir,'BR-CAX_I1PTCLM50_r270v3_sms_amb.clm2.h1.2001-01-01-00000.nc'];...
[dir,'BR-CAX_I1PTCLM50_r270v3_sms_tfe_60.clm2.h1.2001-01-01-00000.nc']...
};


b = ncread(['/Users/kennedy/Desktop/james/data/aug25/',...
    'BR-CAX_I1PTCLM50_r251off_k5g7.clm2.h1.2001-01-01-00000.nc']...
    ,'H2OSOI');
c = zeros(20,52561);
c(:,:)=b(1,:,:);
pore = max(c,[],2);

%dir = '../data/aug25/';
%files = {...
%    [dir,'BR-CAX_I1PTCLM50_r251_k5g7.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r251tfe_k5g7.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r251off_k5g7.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r251offtfe_k5g7.clm2.h1.2001-01-01-00000.nc']};

a=ncinfo(files{1});

offset = 10;
ns     = 20;
nx     = length(files);

if ~exist('fctr','var')
    a          = getvars( files{1} , offset, ns, 1);
    varlist    = {'FCTR','FPSN','BTRAN','VEGWP','SMP','QROOTSINK',...
        ...              1      2      3       4        5     6
        'FCEV','FSH','KSR','GSSUN','GSSHA','ELAI','FGEV','H2OSOI'};
    %     7      8      9     10      11     12     13      14
    vard       = ones(length(varlist),length(files));
    vard(3,:)  = [0,0,1,1];
    vard(4,:)  = [4,4,0,0];
    vard(5,:)  = ns;
    vard(6,:)  = ns;
    vard(9,:)  = [ns,ns,0,0];
    vard(14,:) = ns;
    x          = getmore( files,offset,n,varlist,vard );
end

%targ  = fctr+fgev+fcev;
%vvals = [0:0.25:2,max(vpd)];
tstr={'phs-on','phs-on TFE','phs-off','phs-off TFE'};

oneto = 1:n;

zs=[     0  , 0.0200  , 0.0600  , 0.1200  , 0.2000  , 0.3200  , 0.4800,...
    0.6800  , 0.9200  , 1.2000  , 1.5200  , 1.8800  , 2.2800  , 2.7200,...
    3.2600  , 3.9000  , 4.6400  , 5.4800  , 6.4200  , 7.4600  , 8.6000];
z = zs(1:20)+zs(2:21);
z = z/2;

zr=[1.40e-2,2.73e-2,3.96e-2,5.02e-2,7.02e-2,...
    8.49e-2,9.36e-2,9.62e-2,9.36e-2,8.67e-2,...
    7.68e-2,6.54e-2,5.36e-2,4.67e-2,3.67e-2,...
    2.62e-2,1.71e-2,1.03e-2,5.70e-3,2.92e-3];


f = '../data/daily_sapflow_control.txt';
a=csvread(f,1,0);
p = a(:,2);
d = a(:,1);
f = '../data/daily_sapflow_drought.txt';
a=csvread(f,1,0);
p = [p,a(:,2)];



%************************************************************************
%------------------------------------------------------------------------

ff = [0,2,0,0,0,...
    0,0,0,0,0,...
    0,0,0,0,0,...
    0,0,0,0,0,...
    0];

%2  = water potential
%3  = timeseries
%4  = diurnal
%5  = stress vs. vpd
%6  = conductances
%7  = SON, daytime, qrootsink
%8  = FMA, daytime, qrootsink
%9  = total HR
%10 = nighttime sink
%11 = smp profile
%12 = h2osoi profile
%13 = top3m
%16 = retention curves
%17 = nov03, lwp
%18 = new conductance angle 1
%19 = new conductance angle 2
%20 = look at soil profiles in time
%21 = ibid

if ff(13)>0
xdk = figure;
for i=1:4
    zv = zs(2:end)-zs(1:end-1);
    x=zv(1:13)*h2osoi((1:13)+(i-1)*20,:)+h2osoi(14+(i-1)*20,:)*(3-zs(14));
    plot(1000*x,'LineWidth',1.5)
    hold on
end
legend({'phs-amb','phs-tfe','sms-amb','sms-tfe'},'location','eastoutside')
xlim([0,3*48*365])
set(gca,'xtick',0:48*365:3*48*365)
set(gca,'xticklabel',2001:2004)
xlabel('Time')
title('Water content of Top 3m of Soil Column')
ylabel('mm H2O')

xdk.Units = 'inches';
xdk.Position = [2,2,6,3];
xdk.PaperSize = [6,3];
xdk.PaperPosition = [0,0,6,3];

print(xdk,'../figs2/top3m','-dpdf')
end


if ff(2)>0
    
    %figure 2, SON2013 diurnal mean of veg water potential
    
    % calculate monthly mean midday lwp
    tt=25;
    ixd = mcsec >diurn(tt)&mcsec<diurn(tt+4);
    
    for i=1:2
    lwp(i,:) = 1/101972*splitapply(@mean,vegwp(1+(i-1)*4,ixd),month(ixd)+(year(ixd)-2001)*12);
    end
    

    
    % calculate SON-2003 diurnal mean
    ix  = year==2003&month>8&month<12;
    g   = findgroups(mcsec);
    out = zeros(8,48);
    for i=1:8
        out(i,:) = splitapply(@mean,vegwp(i,ix),g(ix));
    end
    
    %plotting
    xdk = figure;
    x = 0.25:0.5:24;
    xf = 1/101972;  %converts mm to MPa
    
    subplot('position',[0.1, 0.59, 0.43, 0.36])
    plot(x,xf*out(1:4,:)')
    xlim([0 24])
    ylim([-3 0])
    xlabel('Hour')
    ylabel({'Water Potential';'(MPa)'})
    title('AMB')
    text(1,-2.5,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:6:24)
    
    subplot('position',[0.55, 0.59, 0.43, 0.36])
    plot(x,xf*out(5:8,:)')
    set(gca,'xtick',0:6:24)
    ylim([-3 0])
    xlim([0 24])
    xlabel('Hour')
    title('TFE')
    set(gca,'yticklabel',[])
    l = legend('sun-leaf','shade-leaf','stem','root','location','southeast');
    l.Position(1) = 0.45;
    l.Position(2) = 0.65;
    set(gca,'xtick',0:6:24)
    text(21,-2.5,'(b)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.1, 0.1, 0.88, 0.36])
    plot(lwp(1,:),'k-','LineWidth',2)
    hold on
    plot(lwp(2,:),'k:','LineWidth',2)
    ylim([-3,-1])
    xlim([0,36])
    set(gca,'ytick',-3:0.5:-1)
    set(gca,'yticklabel',{'-3','','-2','','-1'})
        set(gca,'xtick',3:3:36)
    %set(gca,'xticklabel',repmat(3:3:12,1,3))
    xlabel('Month')
    ylabel({'Midday Leaf Water Potential';'(MPa)'})
    legend('AMB','TFE','Location','Southwest')
    text(33.8,-1.3,'(c)','FontSize',14,'FontWeight','bold')
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(2)>1
        print(xdk,'../figs2/fig2','-dpdf')
    end
    
    %predawn water potentials?
    out         = xf*out;
    
    disp('AMB MIDDAY')
    disp(out(1:3,27)')
    
    disp('TFE MIDDAY')
    disp(out(5:7,27)')
    
    disp('AMB PREDAWN ROOT')
    disp(out(4,11))
    
    disp('TFE PREDAWN ROOT')
    disp(out(8,11))
    
    disp('delta')
    disp(out(4,11)-out(8,11))
    
    disp('AMB MD ROOT & drop')
    disp([out(4,27),out(4,27)-out(4,11)])
    
    disp('TFE MD ROOT & drop')
    disp([out(8,27),out(8,27)-out(8,11)])
    
    disp('delta drop')
    disp(out(4,27)-out(4,11)-(out(8,27)-out(8,11)))
    
    disp('AMB MD SUN & drop')
    disp([out(1,27),out(4,27)-out(1,27)])
    
    disp('TFE MD SUN & drop')
    disp([out(5,27),out(8,27)-out(5,27)])
    
    disp('delta drop')
    disp(out(4,27)-out(1,27)-(out(8,27)-out(5,27)))
    
    disp('net leaf drop')
    disp(out(1,27)-out(5,27))
    
end


if ff(3)>0
    out = nan(14,36);
    
    % look at timeseries of btran
    % for vegwp, average 12-2 (i.e. 5 timesteps)
    ix  = mcsec>=diurn(25)&mcsec<diurn(29);
    targ = 2.^-((vegwp(1,:)/-250000).^3.95);
    out(1,:) = splitapply(@mean,targ(ix),month(ix)+(year(ix)-2001)*12);
    ix  = mcsec>=diurn(25)&mcsec<diurn(29);
    targ = 2.^-((vegwp(5,:)/-250000).^3.95);
    out(2,:) = splitapply(@mean,targ(ix),month(ix)+(year(ix)-2001)*12);
    
    % for btran, straight average
    out(3,:) = splitapply(@mean,btran(1,:),month+(year-2001)*12);
    out(4,:) = splitapply(@mean,btran(2,:),month+(year-2001)*12);
    
    % GPP
    gppxf = 86400*12/1e6; %umol/m2/s->g/m2/d
    for i=1:4
        out(i+4,:) = gppxf*splitapply(@mean,fpsn(i,:),month+(year-2001)*12);
    end
    
    % transpiration
    for i=1:4
        out(i+8,:) = splitapply(@mean,fctr(i,:),month+(year-2001)*12);
    end
    
    %sapflux
    oneto = (1:730)';
    mvals = [0,cumsum(repmat(eomday(2001,1:12),1,2))];
    out2 = nan(2,24);
    for i=1:2
        x = nan*p(:,1);
        ix = p(:,i)>0;
        x(ix) = p(ix,i);
        for mm=1:24
            ix = oneto>mvals(mm)&oneto<mvals(mm+1);
            if sum(~isnan(x(ix)))>4
                out2(i,mm) = nanmean(x(ix));
            end
        end
    end
    out2 = out2/86400/4e-7; %convert mm/d to W/m2
    out(13:14,13:36) = out2;
    
    %plotting
    xdk = figure;
    
    subplot('Position',[0.07,0.69,0.44,0.28])
    hold on
    plot(out(1,:),'k-','LineWidth',1.5)
    plot(out(2,:),'k:','LineWidth',1.5)
    ylim([0 1])
    xlim([0 36])
    grid on
    ylabel('Stress function')
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    text(2,0.167,'(a)','FontSize',14,'FontWeight','bold')
    title('PHS')
    
    subplot('Position',[0.535,0.69,0.44,0.28])
    hold on
    plot(out(3,:),'k-','LineWidth',1.5)
    plot(out(4,:),'k:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    ylim([0 1])
    xlim([0 36])
    grid on
    text(2,0.167,'(b)','FontSize',14,'FontWeight','bold')
    title('SMS')
    
    subplot('Position',[0.07,0.38,0.44,0.28])
    hold on
    plot(out(5,:),'k-','LineWidth',1.5)
    plot(out(6,:),'k:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    ylim([0 10])
    xlim([0 36])
    grid on
    ylabel('GPP (g/m2/d)')
    text(2,1.67,'(c)','FontSize',14,'FontWeight','bold')
    legend({'AMB','TFE'},'location','SouthEast')
    
    subplot('Position',[0.535,0.38,0.44,0.28])
    hold on
    plot(out(7,:),'k-','LineWidth',1.5)
    plot(out(8,:),'k:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    ylim([0 10])
    xlim([0 36])
    grid on
    text(2,1.67,'(d)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.07,0.07,0.44,0.28])
    hold on
    plot(out(9,:),'k-','LineWidth',1.5)
    plot(out(10,:),'k:','LineWidth',1.5)
    plot(out(13,:),'r-','LineWidth',1.5)
    plot(out(14,:),'r:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    %set(gca,'xticklabel',repmat(3:3:12,1,3))
    ylim([0 150])
    xlim([0 36])
    grid on
    ylabel('T (W/m2)')
    xlabel('Month')
    text(2,25,'(e)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.535,0.07,0.44,0.28])
    hold on
    plot(out(11,:),'k-','LineWidth',1.5)
    plot(out(12,:),'k:','LineWidth',1.5)
    plot(out(13,:),'r-','LineWidth',1.5)
    plot(out(14,:),'r:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    %set(gca,'xticklabel',repmat(3:3:12,1,3))
    set(gca,'yticklabel',[])
    ylim([0 150])
    xlim([0 36])
    grid on
    xlabel('Month')
    text(2,25,'(f)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,6];
    xdk.PaperSize = [7,6];
    xdk.PaperPosition = [0,0,7,6];
    
    if ff(3)>1
        print(xdk,'../figs2/fig3','-dpdf')
    end
    
    ix = year==2003;
    gpp03 = 1800*12*sum(fpsn(:,ix),2)*10^-6;
    gpp03(1)-gpp03(3)
    t03   = 1800*sum(fctr(:,ix),2)*4e-7
    t03(1)-t03(3)
    
end


if ff(4)>0
    
    out = zeros(48,4);
    
    targ = 0*smp(1:4,:);
    targ(1,:) = 2.^-((vegwp(1,:)/-250000).^3.95);
    targ(2,:) = 2.^-((vegwp(5,:)/-250000).^3.95);
    targ(3,:) = btran(1,:);
    targ(4,:) = btran(2,:);
    
    g = findgroups(mcsec);
    for ee=1:4
        ix = (month==9|month==10|month==11) & year==2003;
        %       if ee==2||ee==4
        %           ix = (month==9|month==10|month==11) & year>2001;
        %       else
        %           ix = (month==9|month==10|month==11);
        %       end
        out(:,ee) = splitapply(@mean,targ(ee,ix),g(ix));
    end
    
    xdk = figure;
    
    subplot('Position',[0.08,0.11,0.45,0.84])
    hold on
    plot(0.25:0.5:24,out(:,1),'k-','LineWidth',2)
    plot(0.25:0.5:24,out(:,2),'k:','LineWidth',2)
    set(gca,'ytick',0:0.25:1)
    set(gca,'xtick',0:6:24)
    xlabel('Hour')
    ylabel('Stress function')
    legend('AMB','TFE','Location','SouthEast')
    text(1.2,0.1,'(a)','FontSize',14,'FontWeight','bold')
    ylim([0 1])
    
    subplot('Position',[0.54,0.11,0.45,0.84])
    hold on
    plot(0.25:0.5:24,out(:,3),'k-','LineWidth',2)
    plot(0.25:0.5:24,out(:,4),'k:','LineWidth',2)
    set(gca,'xtick',0:6:24)
    set(gca,'ytick',0:0.25:1)
    set(gca,'yticklabel',[])
    xlabel('Hour')
    text(22,0.1,'(b)','FontSize',14,'FontWeight','bold')
    ylim([0 1])
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    
    if ff(4)>1
        print(xdk,'../figs2/fig4','-dpdf')
    end
    
    % don't time average!
    ix = (month==9|month==10|month==11) & year==2003;
    x = mean(targ(:,ix),2);
    gppxf = 86400*12*1e-6; %umol/m2/s -> g/m2/day
    y = gppxf*mean(fpsn(:,ix),2);
    e = mean(fctr(:,ix)+fgev(:,ix)+fcev(:,ix),2);
    z = targ(:,ix)*fpsn(3,ix)'/sum(fpsn(3,ix));

end


if ff(5)>0
    
    ix1 = smp(1:4,:)==1;
    ix2 = smp(1:4,:)==1;
    ix3 = smp(1:4,:)==1;
    
    x    = vpd;
    y    = 0*smp(1:4,:);
    ffix = fsds>400&fsds<425;
    
    %phs-on, amb
    ee      = 1;
    y(ee,:) = 2.^-((vegwp(1,:)/-250000).^3.95);
    tmp     = repmat(vegwp(4,mcsec==diurn(10)),48,1);
    z       = tmp(1:length(fsds));
    q       = quantile(z(ffix),[1/3,2/3]);
    
    ix1(ee,:) = z<q(1)&ffix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ffix;
    ix3(ee,:) = z>=q(2)&ffix;
    
    subplot(1,2,1)
    plot(z(ffix),y(1,ffix),'.')
    
    %phs-on, tfe
    ee      = 2;
    y(ee,:) = 2.^-((vegwp(5,:)/-250000).^3.95);
    tmp     = repmat(vegwp(8,mcsec==diurn(10)),48,1);
    z       = tmp(1:length(fsds));
    ix      = ffix&year>2001;
    q       = quantile(z(ix),[1/3,2/3]);
    
    ix1(ee,:) = z<q(1)&ix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ix;
    ix3(ee,:) = z>=q(2)&ix;
    
    subplot(1,2,2)
    plot(z(ffix),y(2,ffix),'.')
    
    
    %phs-off, amb
    ll = 41:60;ee=3;
    y(ee,:) = btran(1,:);
    z = zr*smp(ll,:);
    q = quantile(z(ffix),[1/3,2/3]);
    
    
    ix1(ee,:) = z<q(1)&ffix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ffix;
    ix3(ee,:) = z>=q(2)&ffix;
    
    %phs-off, amb
    ll = 61:80;ee=4;
    y(ee,:) = btran(2,:);
    z = zr*smp(ll,:);
    ix = ffix&year>2001;
    q = quantile(z(ix),[1/3,2/3]);
    
    ix1(ee,:) = z<q(1)&ix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ix;
    ix3(ee,:) = z>=q(2)&ix;
    
    %plotting
    xdk = figure;
    
    
    subplot('position',[0.08, 0.56, 0.43, 0.39])
    ee  = 1;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'xticklabel',[])
    set(gca,'ytick',0:0.25:1)
    ylabel('Stress function')
    text(0.1,0.15,'(a)','FontSize',14,'FontWeight','bold')
    title('PHS')
    
    
    subplot('position',[0.08, 0.12, 0.43, 0.39])
    ee  = 2;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
        set(gca,'ytick',0:0.25:1)
    ylabel('Stress function')
    xlabel('VPD (kPa)')
    text(0.1,0.15,'(c)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    ee  = 3;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'ytick',0:0.25:1)
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
title('SMS')
    text(2.9,0.15,'(b)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    ee  = 4;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'ytick',0:0.25:1)
    set(gca,'yticklabel',[])
    xlabel('VPD (kPa)')
    text(2.9,0.15,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(5)>1
        print(xdk,'../figs2/fig5','-dpdf')
    end
    
end



if ff(6) >0
    xdk = figure;
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>4;
    kon(1:40,ix) = ksr(1:40,ix);
    
    
    for ss = 41:80
        
        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-254900;
        kon(ss,ix1) = 0;
    end
    
    out(1:80,1) = nanmedian(kon,2);
    out(1:80,3) = quantile(kon',0.75)-out(:,1)';
    out(1:80,2) = quantile(kon',0.25)-out(:,1)';
    
    subplot('position',[0.07, 0.56, 0.43, 0.39])
    hold on
    plot(1:20,out(1:20,1),'x')
    errorbar(1:20,out(1:20,1),out(1:20,2),out(1:20,3),'x','Color',[0.7 0.1 0.1],'Marker','none')
    set(gca,'xticklabel',[])
    ylabel('Conductance (s-1)')
    ylim([0 7e-9])
    text(17,6e-9,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    hold on
    ll = 21:40;
    plot(1:20,out(ll,1),'x')
    errorbar(1:20,out(ll,1),out(ll,2),out(ll,3),'x','Color',[0.7 0.1 0.1],'Marker','none')
    set(gca,'xticklabel',[])
    box off
    ylim([0 7e-9])
    text(17,6e-9,'(b)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.07, 0.12, 0.43, 0.39])
    hold on
    ll = 41:60;
    plot(1:20,out(ll,1),'x')
    errorbar(1:20,out(ll,1),out(ll,2),out(ll,3),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlabel('Soil Layer')
    ylabel('Conductance (s-1)')
    ylim([0 9e-11])
    text(17,6/7*9e-11,'(c)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    hold on
    ll = 61:80;
    plot(1:20,out(ll,1),'x')
    errorbar(1:20,out(ll,1),out(ll,2),out(ll,3),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlabel('Soil Layer')
    ylim([0 9e-11])
    text(17,6/7*9e-11,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(6)>1
        print(xdk,'../figs2/fig6','-dpdf')
    end
    
    
    mm=min(month(month>8&year==2003&fctr(1,:)>4&kon(67,:)==0));
    dd=min(day(month==mm&year==2003&fctr(1,:)>4&kon(67,:)==0));
    disp('first zero conductance')
    disp([num2str(mm),'-',num2str(dd)])
    
    figure
    ix = month>8&year==2003&fctr(1,:)>4;
    plot(kon(67,ix))
    
end


%daytime SON qrootsink
if ff(7)>0
    out     = zeros(80,4);
    ix_dry  = (mcsec>diurn(12)&mcsec<diurn(37))&(month==9|month==10|month==11)&year==2003;
    
    for i=1:4
        ll = (1:20)+(i-1)*20;
        oo = (1:20)+(i-1)*20;
        
        out(oo,1) =     mean(qrootsink(ll,ix_dry),2);
        out(oo,2) =   median(qrootsink(ll,ix_dry),2);
        out(oo,3) = quantile(qrootsink(ll,ix_dry)',0.25)-out(oo,2)';
        out(oo,4) = quantile(qrootsink(ll,ix_dry)',0.75)-out(oo,2)';
        
    end
    
    xdk = figure;
    
    subplot('position',[0.1, 0.55, 0.44, 0.39])
    ll = 1:20;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-2.5e-6,20e-6])
    set(gca,'xticklabel',[])
    text(18.5,17e-6,'(a)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel({'PHS Soil sink';'(mm/s)'})
    title('AMB')
    
    subplot('position',[0.545, 0.55, 0.44, 0.39])
    ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-2.5e-6,20e-6])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,17e-6,'(b)','FontSize',14,'FontWeight','bold')
    title('TFE')
    
    subplot('position',[0.1, 0.11, 0.44, 0.39])
    ll = 41:60;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    set(gca,'xticklabel',{'',5:5:20})
    ylim([-2.5e-6,20e-6])
    xlabel('Soil Layer')
    text(18.5,17e-6,'(c)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel({'SMS Soil sink';'(mm/s)'})
    
    
    subplot('position',[0.545, 0.11, 0.44, 0.39])
    ll = 61:80;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-2.5e-6,20e-6])
    xlabel('Soil Layer')
    set(gca,'xticklabel',{'',5:5:20})
    set(gca,'yticklabel',[])
    text(18.5,17e-6,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(7)>1
        print(xdk,'../figs2/fig7','-dpdf')
    end
    
    disp('phs-on layer 3 extraction')
    disp(1800*sum(qrootsink(3,ix_dry)))
    
    disp('phs-on layer 8 extraction')
    disp(1800*sum(qrootsink(48,ix_dry)))
    
    disp('tfe deep extraction')
    disp([1800*sum(sum(qrootsink(35:40,ix_dry))),1800*sum(sum(qrootsink(75:80,ix_dry)))])
    
    [a,b]=max(max(qrootsink(21:40,ix_dry),[],2));
    disp('max phs root sink')
    disp(['occurs from layer',num2str(b)])
    disp(a)
    
    p  = max(max(h2osoi(b:20:80,:)));
    dz = 1000*(zs(b+1)-zs(b));
    disp('phs time')
    disp(p*dz/a/3600)
    
    p  = max(reshape(max(h2osoi,[],2),20,4),[],2);
    a  = max(qrootsink(61:80,ix_dry),[],2);
    [x,b] = max(a);
    
    disp('max off root sink')
    disp(['from layer',num2str(b)])
    disp(x)
    
    disp('time')
    [x,b] = min(z'*1000.*p./a/3600);
    disp(['from layer',num2str(b)])
    disp(x)
    disp(a(b))
 
end

if ff(8)>0
    out     = zeros(80,4);
    ix_wet  = (mcsec>diurn(12)&mcsec<diurn(37))&(month==2|month==3|month==4)&year==2003;
    
    for i=1:4
        ll = (1:20)+(i-1)*20;
        oo = (1:20)+(i-1)*20;
        
        out(oo,1) =     mean(qrootsink(ll,ix_wet),2);
        out(oo,2) =   median(qrootsink(ll,ix_wet),2);
        out(oo,3) = quantile(qrootsink(ll,ix_wet)',0.25)-out(oo,2)';
        out(oo,4) = quantile(qrootsink(ll,ix_wet)',0.75)-out(oo,2)';
        
    end
    
    xdk = figure;
    
    subplot('position',[0.1, 0.55, 0.44, 0.39])
    ll = 1:20;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,4e-5])
    set(gca,'xticklabel',[])
    text(18.5,25e-6,'(a)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel({'PHS soil sink','(mm/s)'})
    title('AMB')
    
    subplot('position',[0.545, 0.55, 0.44, 0.39])
    ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,4e-5])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,25e-6,'(b)','FontSize',14,'FontWeight','bold')
    title('TFE')
    
    
    subplot('position',[0.1, 0.11, 0.44, 0.39])
    ll = 41:60;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    set(gca,'xticklabel',{'',5:5:20})
    ylim([-0.5e-5,4e-5])
    xlabel('Soil Layer')
    text(18.5,25e-6,'(c)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel({'SMS soil sink';'(mm/s)'})
    
    
    subplot('position',[0.545, 0.11, 0.44, 0.39])
    ll = 61:80;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,4e-5])
    xlabel('Soil Layer')
    set(gca,'xticklabel',{'',5:5:20})
    set(gca,'yticklabel',[])
    text(18.5,25e-6,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(8)>1
        print(xdk,'../figs2/fig8','-dpdf')
    end
    
    % max tfe qrootsink
    [a,b]=max(max(qrootsink(21:40,ix_wet),[],2));
    disp('max tfe qrootsink')
    disp(['from layer: ',num2str(b)])
    disp(a)
    
    
    % total transpiration
    ix = year==2003&month>1&month<5;
    disp('total tfe transpiration on/off')
    disp([1800*sum(sum(qrootsink(21:40,ix))),...
    1800*sum(sum(qrootsink(61:80,ix)))])
    
    
    
end




if ff(9)>0
    
    hr = qrootsink(1:40,:);
    for ss=1:40
        hr(ss,:) = -qrootsink(ss,:).*(qrootsink(ss,:)<0);
    end
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003;
    ix2 = (~(mcsec<diurn(13)|mcsec>diurn(36)))&year==2003;
    hr1 = [1800*splitapply(@sum,sum(hr(1:20,ix)),month(ix))',...
        1800*splitapply(@sum,sum(hr(1:20,ix2)),month(ix2))'];
    hr2 = [1800*splitapply(@sum,sum(hr(21:40,ix)),month(ix))',...
        1800*splitapply(@sum,sum(hr(21:40,ix2)),month(ix2))'];
    
    
    xdk = figure;
    
    subplot('Position',[0.08,0.11,0.45,0.84])
    b=bar(hr1,'stacked');
    b(1).FaceColor = [0.4,0.4,0.6];
    b(2).FaceColor = [0.6,0.6,0.8];
    b(1).EdgeColor = [0.35,0.35,0.6];
    b(2).EdgeColor = [0.55,0.55,0.8];
    xlabel('Month')
    ylabel('Total HR (mm)')
    title('AMB')
    xlim([0 13])
    ylim([0 80])
    text(1,90,'(a)','FontSize',14,'FontWeight','bold')
    
    
    subplot('Position',[0.54,0.11,0.45,0.84])
    b=bar(hr2,'stacked');
    b(1).FaceColor = [0.4,0.4,0.6];
    b(2).FaceColor = [0.6,0.6,0.8];
    b(1).EdgeColor = [0.35,0.35,0.6];
    b(2).EdgeColor = [0.55,0.55,0.8];
    %bar(hr2,'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    xlabel('Month')
    xlim([0 13])
    ylim([0 80])
    title('TFE')
    set(gca,'yticklabel',[])
    text(11,90,'(b)','FontSize',14,'FontWeight','bold')
    legend({'Nighttime','Daytime'})
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(9)>1
        print(xdk,'../figs2/fig9','-dpdf')
    end
    
    disp('Total HR')
    disp([sum(sum(hr1)),sum(sum(hr2))])
    
    disp('Day/night')
    disp([sum(hr1),sum(hr2)])
    
    disp('% day')
    disp([sum(hr1(:,2))/sum(sum(hr1)),sum(hr2(:,2))/sum(sum(hr2))])
    
    disp('corrs')
    p = splitapply(@sum,1800*prec(year==2003),month(year==2003))';
    disp([corr(sum(hr1,2),p),corr(sum(hr2,2),p)])
    
    disp('amb higher these months')
    disp([(1:12)',sum(hr1,2)>sum(hr2,2)])
    
    ix  = year==2003&mcsec==diurn(11);
    x   = splitapply(@mean,vegwp(4,ix),month(ix));
    disp('amb corr with pd vwp')
    disp(corr(x',sum(hr1,2)))
    lm0 = fitlm(p,sum(hr1,2));
    lm1 = fitlm([p,x'],sum(hr1,2));
    x   = splitapply(@mean,vegwp(8,ix),month(ix));
    disp('amb corr with pd vwp')
    disp(corr(x',sum(hr2,2)))
    lm2 = fitlm([p,x'],sum(hr2,2));
    lm3 = fitlm(p,sum(hr2,2));
end


if ff(10)>0
    
    out = zeros(80,4);
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003&(month==9|month==10|month==11);  %night time
    ll = 1:20;
    oo = 1:20;
    out(oo,1) =     mean(qrootsink(ll,ix),2);
    out(oo,2) =   median(qrootsink(ll,ix),2);
    out(oo,3) = quantile(qrootsink(ll,ix)',0.25)-out(oo,2)';
    out(oo,4) = quantile(qrootsink(ll,ix)',0.75)-out(oo,2)';
    
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003&(month==2|month==3|month==4);  %night time
    oo = 21:40;
    out(oo,1) =     mean(qrootsink(ll,ix),2);
    out(oo,2) =   median(qrootsink(ll,ix),2);
    out(oo,3) = quantile(qrootsink(ll,ix)',0.25)-out(oo,2)';
    out(oo,4) = quantile(qrootsink(ll,ix)',0.75)-out(oo,2)';
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003&(month==9|month==10|month==11);  %night time
    ll = 21:40;
    oo = 41:60;
    out(oo,1) =     mean(qrootsink(ll,ix),2);
    out(oo,2) =   median(qrootsink(ll,ix),2);
    out(oo,3) = quantile(qrootsink(ll,ix)',0.25)-out(oo,2)';
    out(oo,4) = quantile(qrootsink(ll,ix)',0.75)-out(oo,2)';
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003&(month==2|month==3|month==4);  %night time
    ll = 21:40;
    oo = 61:80;
    out(oo,1) =     mean(qrootsink(ll,ix),2);
    out(oo,2) =   median(qrootsink(ll,ix),2);
    out(oo,3) = quantile(qrootsink(ll,ix)',0.25)-out(oo,2)';
    out(oo,4) = quantile(qrootsink(ll,ix)',0.75)-out(oo,2)';
    
    xdk = figure;
    
    subplot('position',[0.08, 0.55, 0.45, 0.39])
    ll = 1:20;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    set(gca,'xticklabel',[])
    set(gca,'ytick',-0.4e-5:0.4e-5:2e-5)
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel('Soil sink (mm/s)')
    ylim([-0.6e-5,0.81e-5])
    text(18.5,0.66e-5,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.55, 0.45, 0.39])
    ll = 41:60;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.6e-5,0.81e-5])
    set(gca,'ytick',-0.4e-5:0.4e-5:2e-5)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,0.66e-5,'(b)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.08, 0.11, 0.45, 0.39])
    
    ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    
    xlabel('Soil Layer')
    ylabel('Soil sink (mm/s)')
    
    ylim([-1e-5,2.5e-5])
    text(18.5,2.1e-5,'(c)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.11, 0.45, 0.39])
    ll = 61:80;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-1e-5,2.5e-5])
    xlabel('Soil Layer')
    set(gca,'yticklabel',[])
    text(18.5,2.1e-5,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    
    if ff(10)>1
        print(xdk,'../figs2/fig10','-dpdf')
    end
    
    ix_dry = year==2003&month>8&month<12&(mcsec<diurn(13)|mcsec>diurn(36));
    
    x = sum(qrootsink(1:20,ix_dry),2);
    a = 1800*sum(x(x<0))
    1800*sum(x(x>0))
    
    x = out(1:20,1);
    disp('AMB dry net HR')
    disp(1800*24*91*sum(x(x<0)))
    
    disp('TFE dry net HR')
    x = out(41:60,1);
    disp(1800*24*91*sum(x(x<0)))
    
    disp('Amb dry deep extraction')
    disp(1800*24*91*sum(out(15:20,1)))
    
    disp('TFE dry deep extraction')
    disp(1800*24*91*sum(out(55:60,1)))
    
    x = out(21:40,1);
    disp('AMB wet net HR')
    disp(1800*24*91*sum(x(x<0)))
    
    disp('TFE wet net HR')
    x = out(61:80,1);
    disp(1800*24*91*sum(x(x<0)))
    
   
end


if ff(11)>0
    
    x=mean(smp(:,year==2003&month>8&month<12),2)/101972;
    y=mean(smp(:,year==2003&month>1&month<5),2)/101972;
    
    
    out = zeros(40,8);
    for i=1:4
    tmp = repmat(x((1:20)+(i-1)*20)',2,1);
    out(:,i) = tmp(:);
    end
    for i=1:4
    tmp = repmat(y((1:20)+(i-1)*20)',2,1);
    out(:,i+4) = tmp(:);
    end
    
    zv = zeros(40,1);
    zv(1) = 0;
    ct = 1;
    for i = 1:19
        ct = ct+1;
        zv(ct) = zs(i+1)-.0001;
        ct = ct+1;
        zv(ct) = zs(i+1)+.0001;
    end
    zv(end) = zs(end);
    
    xdk = figure;
    
    subplot('Position',[0.06,0.56,0.42,0.42])
    plot(out(:,5),-zv,'k','LineWidth',1.5)
    hold on
    plot(out(:,6),-zv,'-','LineWidth',1.5,'Color',[0.7,0.7,0.7])
    xlim([-2.5,0])
    ylim([-8.5,0])
    text(-2.35,-7.3,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xticklabel',[])
    set(gca,'ytick',-8:2:0)
    set(gca,'yticklabel',8:-2:0)
    box off
    ylabel('Depth (m)')
    yyaxis right
    ylim([-10,0])
    set(gca,'ytick',-z(20:-5:5))
    set(gca,'yticklabel',[])
    ax = gca;
    ax.YColor = [0,0,0];
    legend('AMB','TFE','Location','NorthWest')
    
    subplot('Position',[0.51,0.56,0.42,0.42])
    plot(out(:,1),-zv,'k','LineWidth',1.5)
    hold on
    plot(out(:,2),-zv,'-','LineWidth',1.5,'Color',[0.7,0.7,0.7])
    xlim([-2.5,0])
    ylim([-8.5,0])
    text(-2.35,-7.3,'(b)','FontSize',14,'FontWeight','bold')
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    box off
    yyaxis right
    ylim([-10,0])
    set(gca,'ytick',-z(20:-5:5))
    set(gca,'yticklabel',[])
    ax = gca;
    ax.YColor = [0,0,0];
    ylabel('Soil Layer')
    set(gca,'yticklabel',20:-5:5)
    
    subplot('Position',[0.06,0.11,0.42,0.42])
    plot(out(:,7),-zv,'k','LineWidth',1.5)
    hold on
    plot(out(:,8),-zv,'-','LineWidth',1.5,'Color',[0.7,0.7,0.7])
    xlim([-2.5,0])
    ylim([-8.5,0])
    text(-2.35,-7.3,'(c)','FontSize',14,'FontWeight','bold')
    set(gca,'ytick',-8:2:0)
    set(gca,'yticklabel',8:-2:0)
    box off
    ylabel('Depth (m)')
    yyaxis right
    ylim([-10,0])
    set(gca,'ytick',-z(20:-5:5))
    set(gca,'yticklabel',[])
    ax = gca;
    ax.YColor = [0,0,0];
    xlabel('Water Potential (MPa)')
    
    
    subplot('Position',[0.51,0.11,0.42,0.42])
    plot(out(:,3),-zv,'k','LineWidth',1.5)
    hold on
    plot(out(:,4),-zv,'-','LineWidth',1.5,'Color',[0.7,0.7,0.7])
    xlim([-2.5,0])
    ylim([-8.5,0])
    text(-2.35,-7.3,'(d)','FontSize',14,'FontWeight','bold')
    set(gca,'yticklabel',[])
    box off
    yyaxis right
    ylim([-10,0])
    set(gca,'ytick',-z(20:-5:5))
    set(gca,'yticklabel',20:-5:5)
    ax = gca;
    ax.YColor = [0,0,0];
    ylabel('Soil Layer')
    xlabel('Water Potential (MPa)')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(11)>1
        print(xdk,'../figs2/fig11','-dpdf')
    end
    
    ix = year==2003&month>8&month<12;
    disp('PHS off (amb/tfe): ')
    disp([mean(zr*smp(41:60,ix))/101972,...
        mean(zr*smp(61:80,ix))/101972])
    
    ix = year==2003&month>8&month<12&mcsec==diurn(11);
    disp('PHS on (amb/tfe): ')
    disp([mean(vegwp(4,ix))/101972,...
    mean(vegwp(8,ix))/101972])

    targ = fctr+fcev+fgev;
    disp('total ET: ')
    y = 1800*sum(targ,2)*4e-10;
    disp(y)
    
    zv = zs(2:end)-zs(1:end-1);
    x=zv*reshape(h2osoi(:,end),20,4)
    y(3)-y(1)
    x(1)-x(3)
    y(4)-y(2)
    x(2)-x(4)
end

if ff(12)>0
    x=mean(h2osoi(:,year==2003&month>8&month<12),2);
    
    out = zeros(40,4);
    for i=1:4
    tmp = repmat(x((1:20)+(i-1)*20)',2,1);
    out(:,i) = tmp(:);
    end
    
        zv = zeros(40,1);
    zv(1) = 0;
    ct = 1;
    for i = 1:19
        ct = ct+1;
        zv(ct) = zs(i+1)-.0001;
        ct = ct+1;
        zv(ct) = zs(i+1)+.0001;
    end
    zv(end) = zs(end);

    subplot(1,2,1)
    plot(out(:,1),-zv,'LineWidth',2)
    hold on
    plot(out(:,3),-zv,'LineWidth',2)
    xlim([0 1])
    title('AMB')
    subplot(1,2,2)
    plot(out(:,2),-zv,'LineWidth',2)
    hold on
    plot(out(:,4),-zv,'LineWidth',2)
    title('TFE')
    xlim([0 1])
    
 
end

if ff(14)>0
    
    yy  = 2002;
    dd  = 200;
    for i = 1:529
        dd = dd+1;
        if dd==366
            dd = 1;
            yy = 2003;
        end
        
        subplot(1,2,1)
        x = mean(smp(61:80,year==yy&doy==dd),2)/101972;
        %x = mean(h2osoi(61:80,year==yy&doy==dd),2);
        plot(x,-1:-1:-20,'k-','LineWidth',2)
        xlim([-2.7,0])
        title(dd)
        subplot(1,2,2)
        x = mean(smp(21:40,year==yy&doy==dd),2)/101972;
        %x = mean(h2osoi(21:40,year==yy&doy==dd),2);
        plot(x,-1:-1:-20,'k-','LineWidth',2)
        xlim([-2.7,0])
        pause(0.1)
    end   
    
    
end



if ff(15)>0
    plot(vegwp(1,mcsec==diurn(10))/101972)
    ylim([-1,0])
    figure
    plot(vegwp(5,mcsec==diurn(10))/101972)
    ylim([-1,0])
end

if ff(16)>0
    s = smp(1:20,24:48*5:end);
    s = s(:)/101.900;
    h = h2osoi(1:20,24:48*5:end);
    h = h(:);
    
    semilogx(-s,h,'.')
    set(gca,'xtick',10.^(0:3))
    set(gca,'xticklabel',-[1,10,100,1000])
    grid on
    hold on
    semilogx(1000,0.15,'rx')
    semilogx(1,0.35,'rx')
end

if ff(17)>0
    ix  = year==2003&month==11;
    out = splitapply(@mean,vegwp(5,ix),findgroups(mcsec(ix)));
    plot(out/101972)
    hold on
    out = splitapply(@mean,vegwp(8,ix),findgroups(mcsec(ix)));
    plot(out/101972)

end


if ff(18)>0
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>4;
    kon(1:40,:) = ksr(1:40,:);
    
    
    for ss = 41:80
        
        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-254900;
        kon(ss,ix1) = 0;
    end
    
    ix = year==2003;
    pday2003 = 1800*splitapply(@sum,prec(ix),doy(ix));
    
    xdk = figure;
     
    subplot('Position',[0.09,0.69,0.9,0.26])
    plot((1:sum(year==2003))/48,kon(3,year==2003),'.','Color',[0.4,0.4,0.4])
    ylabel({'PHS Conductance';'(1/s)'})
    xlim([-1,366])
    set(gca,'xticklabel',[])
    text(5,1e-9,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.09,0.39,0.9,0.26])
    plot((1:sum(year==2003))/48,kon(43,year==2003),'.','Color',[0.4,0.4,0.4])
    ylabel({'SMS Conductance';'(1/s)'})
    xlim([-1,366])
    set(gca,'xticklabel',[])
    text(5,5e-11,'(b)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.09,0.09,0.9,0.26])
    bar(pday2003,'FaceColor',[0.4,0.4,0.4])
    ylabel('Rain (mm/d)')
    xlim([-1,366])
    xlabel('Day of 2003')
    text(5,5/6*80,'(c)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    if ff(18)>1
    print(xdk,'../figs2/fig6','-dpdf')
    end

    
end

if ff(19)>0
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>4;
    kon(1:40,:) = ksr(1:40,:);
    
    
    for ss = 41:80
        
        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-254900;
        kon(ss,ix1) = 0;
    end
    
    ix = year==2003;
    out = zeros(80,365);
    for ss=1:80
    out(ss,:) = splitapply(@nanmean,kon(ss,ix),doy(ix));
    end
    
    ss = 3;
    xdk = figure;
    subplot('Position',[0.1,0.56,0.88,0.39])
    plot(out(ss,:),'k','LineWidth',2)
    hold on
    plot(out(ss+20,:),'k:','LineWidth',2)
    ylabel({'PHS Daily Mean';'Conductance (s^{-1})'})
    xlim([-1,366])
    set(gca,'xticklabel',[])
    text(5,1e-9,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.1,0.11,0.88,0.39])
    plot(out(ss+40,:),'k','LineWidth',2)
    hold on
    plot(out(ss+60,:),'k:','LineWidth',2)
    xlim([-1,366])
    ylabel({'SMS Daily Mean';'Conductance (s^{-1})'})
    legend({'AMB','TFE'},'location','northwest')
    xlabel('Day of 2003')
    text(345,5/6*8e-11,'(b)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches'
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(19)>1
    print(xdk,'../figs2/fig6a','-dpdf')
    end
    
end


if ff(20)>0
    
    yy=2002;
    xdk = figure;
    %resample to height
    x = 0;
    ymax = [-0.8,-2.5];
    vsn  = {'PHS','SMS'};
    pnl  = {'(a)','(b)'};
    for ss=[20,60]
        x = x+1;
    nz  = round(100*(zs(2:end)-zs(1:end-1)));
    out = zeros(860,n);
    ct  = 0; 
    
    for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) = repmat(smp(ss+i,:),nz(i),1)/101972;
            ct = ct+nz(i);
    end
    addpath('/Users/kennedy/Documents/MATLAB/othercolor')
    colormap(othercolor('BrBG4'))
    subplot(2,1,x)
    imagesc((1:n)/48,(1:860)/100,out,[ymax(x),0])
    title(vsn{x})
    ylim([0,8.605])
    ylabel('Depth (m)')
    if x==2
        xlabel('Days since Jan 1, 2001')
    end
    set(gca,'ytick',0:2:8)
    c = colorbar;
    ylabel(c,'Soil Potential (MPa)','FontSize',11)
    text(30,1,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    print(xdk,'../figs2/fig11','-dpdf')
 
    
end

if ff(21)>0

    for dd=1:1
        subplot(2,2,1)
        plot(mean(smp(22:40,year==2003&doy==dd),2)/101972,-2:-1:-20)
        xlim([-1,0])
        subplot(2,2,2)
        plot(mean(smp(62:80,year==2003&doy==dd),2)/101972,-2:-1:-20)
        
        xlim([-2.5,0])
        pause(0.1)

    end

end