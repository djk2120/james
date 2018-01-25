close all

dir = '../data/jan24/';

files = {...
    [dir,'BR-CAX_I1PTCLM50_r270.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r270_on_tfe.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r270_off_amb.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r270_off_tfe.clm2.h1.2001-01-01-00000.nc']};
    

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

ff = [0,0,0,0,2,...
    0,0,0,0,0,...
    0,0,0,0,0,...
    0,0,0];

%2  = water potential
%5  = stress vs. vpd
%6  = conductances
%7  = SON, daytime, qrootsink
%9  = total HR

if ff(2)>0

    %figure 2, SON2013 diurnal mean of veg water potential

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
    subplot('position',[0.08, 0.12, 0.43, 0.83])
    plot(x,xf*out(1:4,:)')
    xlim([0 24])
    ylim([-2.5 0])
    xlabel('Hour')
    ylabel('Water Potential (MPa)')
    title('AMB')
    text(1.5,-2.3,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:6:24)
    subplot('position',[0.54, 0.12, 0.43, 0.83])
    plot(x,xf*out(5:8,:)')
    set(gca,'xtick',0:6:24)
    ylim([-2.5 0])
    xlim([0 24])
    xlabel('Hour')
    title('TFE')
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    set(gca,'yticklabel',[])
    l = legend('sun-leaf','shade-leaf','stem','root','location','southeast');
    l.Position(1) = 0.45;
    l.Position(2) = 0.22;
    set(gca,'xtick',0:6:24)
    text(21,-2.3,'(b)','FontSize',14,'FontWeight','bold')
    
    if ff(2)>1
        print(xdk,'../figs/fig2','-dpdf')
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
    

if ff(18)>0
    
    x=mean(smp(:,year==2003&month>8&month<12),2)/101972;
    y=mean(smp(:,year==2003&month>1&month<5),2)/101972;
    
    xdk = figure;
    
    subplot('Position',[0.07,0.56,0.42,0.42])
    plot(y(1:20),-z,'k','LineWidth',1.5)
    hold on
    plot(y(21:40),-z,'k:','LineWidth',1.5)
    xlim([-2.75,0])
    ylim([-8.5,0])
    text(-2.6,-7.3,'(a)','FontSize',14,'FontWeight','bold')
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
    plot(x(1:20),-z,'k','LineWidth',1.5)
    hold on
    plot(x(21:40),-z,'k:','LineWidth',1.5)
    xlim([-2.75,0])
    ylim([-8.5,0])
    text(-2.6,-7.3,'(b)','FontSize',14,'FontWeight','bold')
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

    subplot('Position',[0.07,0.11,0.42,0.42])
    plot(y(41:60),-z,'k','LineWidth',1.5)
    hold on
    plot(y(61:80),-z,'k:','LineWidth',1.5)
    xlim([-2.75,0])
    ylim([-8.5,0])
    text(-2.6,-7.3,'(c)','FontSize',14,'FontWeight','bold')
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
    plot(x(41:60),-z,'k','LineWidth',1.5)
    hold on
    plot(x(61:80),-z,'k:','LineWidth',1.5)
    xlim([-2.75,0])
    ylim([-8.5,0])
    text(-2.6,-7.3,'(d)','FontSize',14,'FontWeight','bold')
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
    
    if ff(18)>1
        print(xdk,'../figs/fig13','-dpdf')
    end
end


if ff(17)>0
    ix  = mcsec==diurn(10);
    out = splitapply(@mean,vegwp(4,ix),month(ix)+(year(ix)-2001)*12);

        out2 = splitapply(@mean,vegwp(8,ix),month(ix)+(year(ix)-2001)*12);
    plot(out-out2)
    
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
    grid on
    
    ix  = mcsec>=diurn(25)&mcsec<diurn(29);
    out = splitapply(@mean,vegwp(1,ix),month(ix)+(year(ix)-2001)*12);
    figure
    plot(out)
    hold on
    out = splitapply(@mean,vegwp(5,ix),month(ix)+(year(ix)-2001)*12);
    plot(out)
    
    x = -1e4:-1e4:-5e5;
    s = 2.^-((x/-250000).^3.95);
    figure
    plot(x,s)
    grid on
    st1 = 2^-((150000/250000)^3.95)
    st2 = 2^-((160000/250000)^3.95)
    
end

if ff(16)>0
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

    subplot('Position',[0.07,0.69,0.44,0.285])
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
    
    subplot('Position',[0.535,0.69,0.44,0.285])
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
    
    
    subplot('Position',[0.07,0.38,0.44,0.285])
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
    
    subplot('Position',[0.535,0.38,0.44,0.285])
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
    
    subplot('Position',[0.07,0.07,0.44,0.285])
            hold on
    plot(out(9,:),'k-','LineWidth',1.5)
    plot(out(10,:),'k:','LineWidth',1.5)
    plot(out(13,:),'r-','LineWidth',1.5)
    plot(out(14,:),'r:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
    ylim([0 150])
    xlim([0 36])
    grid on
    ylabel('T (W/m2)')
    xlabel('Month')
    text(2,25,'(e)','FontSize',14,'FontWeight','bold')
    
    subplot('Position',[0.535,0.07,0.44,0.285])
            hold on
    plot(out(11,:),'k-','LineWidth',1.5)
    plot(out(12,:),'k:','LineWidth',1.5)
    plot(out(13,:),'r-','LineWidth',1.5)
    plot(out(14,:),'r:','LineWidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
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
    
    if ff(16)>1
        print(xdk,'../figs/fig12','-dpdf')
    end
end


if ff(15)>0
    oneto = (1:730)';
    mvals = [0,cumsum(repmat(eomday(2001,1:12),1,2))];
    out = nan(24,2);
    
    for i=1:2
        x = nan*p(:,1);
        ix = p(:,i)>0;
        x(ix) = p(ix,i);
        
        for mm=1:24
            ix = oneto>mvals(mm)&oneto<mvals(mm+1);
            if sum(~isnan(x(ix)))>4
                out(mm,i) = nanmean(x(ix));
            end
        end
    end
    
    ix    = p(:,1)>0;
    plot(oneto(ix),p(ix,1),'.')
    hold on
    ix    = p(:,2)>0;
    plot(oneto(ix),p(ix,2),'.')
    set(gca,'xtick',cumsum(repmat(eomday(2001,1:12),1,2)))
    grid on
    set(gca,'xticklabel',1:24)
    
    
    mx = 0.5*(mvals(1:end-1)+mvals(2:end));
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(mx,out)
    
end

if ff(14)>0
    
    %calc monthly mean T
    out = zeros(4,36);
    for i=1:4
        out(i,:) = splitapply(@mean,fctr(i,:),month+(year-2001)*12);
    end
    
    %calc monthly mean sapflow
    oneto = (1:730)';
    mvals = [0,cumsum(repmat(eomday(2001,1:12),1,2))];
    out2 = nan(24,2);
    for i=1:2
        x = nan*p(:,1);
        ix = p(:,i)>0;
        x(ix) = p(ix,i);
        for mm=1:24
            ix = oneto>mvals(mm)&oneto<mvals(mm+1);
            if sum(~isnan(x(ix)))>4
                out2(mm,i) = nanmean(x(ix));
            end
        end
    end
    out2 = out2/86400/4e-7; %convert mm/d to W/m2

    xdk = figure;
    
    subplot('position',[0.08, 0.56, 0.43, 0.39])
    hold on
    plot(out(1,:),'k-','Linewidth',1.5)
    plot(out(2,:),'k:','Linewidth',1.5)
    plot(13:36,out2(:,1),'r-','Linewidth',1.5)
    plot(13:36,out2(:,2),'r:','Linewidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    grid on
    xlim([0 36])
    ylim([0 150])
    ylabel('T (W/m2)')
    text(2,25,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    hold on
    plot(out(3,:),'k-','Linewidth',1.5)
    plot(out(4,:),'k:','Linewidth',1.5)
    plot(13:36,out2(:,1),'r-','Linewidth',1.5)
    plot(13:36,out2(:,2),'r:','Linewidth',1.5)
    xlim([0 36])
    ylim([0 150])
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    grid on
    text(2,25,'(b)','FontSize',14,'FontWeight','bold')
    
    out = zeros(4,36);
    for i=1:4
        out(i,:) = splitapply(@mean,fpsn(i,:),month+(year-2001)*12);
    end
    gppxf = 86400*12*1e-6; %umol/m2/s -> g/m2/day
    out   = out*gppxf;
    
    
    subplot('position',[0.08, 0.12, 0.43, 0.39])
    hold on
    plot(out(1,:),'k-','Linewidth',1.5)
    plot(out(2,:),'k:','Linewidth',1.5)
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
    grid on
    xlim([0 36])
    ylim([0 10])
    ylabel('GPP (g/m2/d)')
    text(2,1.67,'(c)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    hold on
    plot(out(3,:),'k-','Linewidth',1.5)
    plot(out(4,:),'k:','Linewidth',1.5)
    xlim([0 36])
    ylim([0 10])
    set(gca,'xtick',3:3:36)
    set(gca,'xticklabel',repmat(3:3:12,1,3))
    set(gca,'yticklabel',[])
    grid on
    text(2,1.67,'(d)','FontSize',14,'FontWeight','bold')
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(14)>1
        print(xdk,'../figs/fig11','-dpdf')
    end
    
end


if ff(13)>0
    ix  = year==2003;
    
    g   = month(ix);
    
    subplot('Position',[0.08 0.56 0.9 0.42])
    out = splitapply(@sum,1800*sum(qrootsink(15:20,ix)),g);
    plot(out,'k-','LineWidth',2)
    hold on
    out = splitapply(@sum,1800*sum(qrootsink(35:40,ix)),g);
    plot(out,'k:','LineWidth',2)
    xlim([0 13])
    ylim([-10 60])
    ylabel('Soil sink (mm)')
    set(gca,'xtick',1:12)
    set(gca,'xticklabel',[])
    
    subplot('Position',[0.08 0.11 0.9 0.42])
    out = splitapply(@sum,1800*sum(qrootsink(55:60,ix)),g);
    plot(out,'k-','LineWidth',2)
    hold on
    out = splitapply(@sum,1800*sum(qrootsink(75:80,ix)),g);
    plot(out,'k:','LineWidth',2)
    xlim([0 13])
    ylim([-10 60])
    xlabel('Month')
    ylabel('Soil sink (mm)')
    set(gca,'xtick',1:12)

    close all
    out = zeros(4,2);
    for i=1:4
        for yy=1:2
        out(i,yy) = 1800*sum(sum(qrootsink((15:20)+(i-1)*20,year==(yy+2001))));
        end
    end
    
    figure
    bar(out)
    
end
    
    

if ff(12)>0
    
   out = zeros(36,4);

   for ee = 1:4
       out(:,ee) = splitapply(@mean,fctr(ee,:),month+(year-2001)*12);
   end
   
   plot(out)
   ylim([0 150])
   legend('onAMB','onTFE','offAMB','offTFE') 
    
    
end

if ff(11)>0
    
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
   
   
   if ff(11)>1
       print(xdk,'../figs/fig10','-dpdf')
   end
   
   % don't time average!
   ix = (month==9|month==10|month==11) & year==2003; 
    x = mean(targ(:,ix),2);
    gppxf = 86400*12*1e-6; %umol/m2/s -> g/m2/day
    y = gppxf*mean(fpsn(:,ix),2);
    
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
    
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    ee  = 2;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9) 
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'ytick',0:0.25:1)
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    text(2.9,0.15,'(b)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.08, 0.12, 0.43, 0.39])
    ee  = 3;
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
        print(xdk,'../figs/fig5','-dpdf')
    end

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
    
    subplot('position',[0.07, 0.55, 0.45, 0.39])
    ll = 1:20;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,3e-5])
    set(gca,'xticklabel',[])
    text(18.5,25e-6,'(a)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel('Soil sink (mm/s)')
    
    subplot('position',[0.535, 0.55, 0.45, 0.39])
    ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,3e-5])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,25e-6,'(b)','FontSize',14,'FontWeight','bold')
    
    
        subplot('position',[0.07, 0.11, 0.45, 0.39])
    ll = 41:60;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    set(gca,'xticklabel',{'',5:5:20})
    ylim([-0.5e-5,3e-5])
    xlabel('Soil Layer')
    text(18.5,25e-6,'(c)','FontSize',14,'FontWeight','bold')
    ax = gca;
    ax.YAxis.Exponent = -5;
    ylabel('Soil sink (mm/s)')
    
    
    subplot('position',[0.535, 0.11, 0.45, 0.39])
    ll = 61:80;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-0.5e-5,3e-5])
    xlabel('Soil Layer')
    set(gca,'xticklabel',{'',5:5:20})
    set(gca,'yticklabel',[])
    text(18.5,25e-6,'(d)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(8)>1
        print(xdk,'../figs/fig8','-dpdf')
    end
    
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
    
    subplot('position',[0.07, 0.55, 0.45, 0.39])
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
    ylabel('Soil sink (mm/s)')
    
    subplot('position',[0.535, 0.55, 0.45, 0.39])
    ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-2.5e-6,20e-6])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,17e-6,'(b)','FontSize',14,'FontWeight','bold')
    
    
        subplot('position',[0.07, 0.11, 0.45, 0.39])
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
    ylabel('Soil sink (mm/s)')
    
    
    subplot('position',[0.535, 0.11, 0.45, 0.39])
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
        print(xdk,'../figs/fig7','-dpdf')
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


if ff(1)>0

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
    
    subplot('position',[0.06, 0.55, 0.45, 0.39])
    ll = 1:20;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    set(gca,'xticklabel',[])
    set(gca,'ytick',-1e-5:1e-5:2e-5)
    ylabel('Soil sink (mm/s)')
    ylim([-1e-5,1e-5])
    text(18.5,0.75e-5,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.53, 0.55, 0.45, 0.39])
    ll = 41:60;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-1e-5,1e-5])
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    text(18.5,0.75e-5,'(b)','FontSize',14,'FontWeight','bold')
    
    
    subplot('position',[0.06, 0.11, 0.45, 0.39])

        ll = 21:40;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    
    xlabel('Soil Layer')
        ylabel('Soil sink (mm/s)')
    
    ylim([-1e-5,3e-5])
    text(18.5,2.5e-5,'(c)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.53, 0.11, 0.45, 0.39])
    ll = 61:80;
    hold on
    bar(out(ll,1),'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    errorbar(1:20,out(ll,2),out(ll,3),out(ll,4),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0 21])
    ylim([-1e-5,3e-5])
    xlabel('Soil Layer')
    set(gca,'yticklabel',[])
    text(18.5,2.5e-5,'(d)','FontSize',14,'FontWeight','bold')

    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];

    
    if ff(1)>1
        print(xdk,'../figs/fig6','-dpdf')
    end
    
    
    x = reshape(out(:,1),20,4);
    x = x(:,[1,3]);
    y = x;
    x = x.*(x<0);
    y = y.*(y>0);
    sum(x)*1800*91*24
    
    figure
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003&(month==9|month==10|month==11);  %night time
    t = doy(ix)+mcsec(ix)/max(diurn);
    x = qrootsink(22,ix);
       
    xv = -1e-5:1e-6:7e-5;
    out = zeros(length(xv),1);
    for i=1:length(xv)-1
        out(i)=sum(x>xv(i)&x<=xv(i+1));        
    end
    plot(xv,out)
    
    
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
    ylim([0 100])
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
    ylim([0 100])
    title('TFE')
    set(gca,'yticklabel',[])
    text(11,90,'(b)','FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(9)>1
        print(xdk,'../figs/fig9','-dpdf')
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
end

if ff(6) >0
    xdk = figure;
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>4;
    kon(1:40,ix) = ksr(1:40,ix);

    
    for ss = 41:80

        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-255000;
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
        print(xdk,'../figs/fig6','-dpdf')
    end
    
    close all
    mm=min(month(year==2003&fctr(1,:)>4&kon(67,:)==0));
    dd=min(day(month==mm&year==2003&fctr(1,:)>4&kon(67,:)==0));
    disp('first zero conductance')
    disp([num2str(mm),'-',num2str(dd)])
    
end



if ff(3)>0
    kon  = nan*smp(1:20,:);
    for i=1:20
        x = vegwp(8,:)-smp(i+20,:)-zs(i+1)*1000;
        ix = fsds>1&abs(x)>1000&year>2001;
        kon(i,ix) = -qrootsink(20+i,ix)./x(ix);
    end
    
    subplot('position',[0.54, 0.12, 0.43, 0.8])
    
    out = nanmedian(kon,2);
    plot(out,'o','Color',[0.6,0.6,0.8])
    x = quantile(kon',[0.25,0.75])';
    x = abs(x-[out,out]);
    hold on
    errorbar(1:20,out,x(1:20,1),x(1:20,2),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0,21])
    xlabel('Soil Layer')
    title('PHS on')
    
    
    box off
    k2 = nan*smp(61:80,:);
    for i=61:80
        
        ix          = year>2001&fsds>1;
        k2(i-60,ix) = qrootsink(i,ix)./min(189000,(smp(i,ix)+255000));
        
        ix1       = year>2001&fsds>1&smp(i,:)<=-255000;
        k2(i-60,ix1) = 0;
        
        
    end
    
    subplot('position',[0.08, 0.12, 0.43, 0.8])
    out = nanmedian(k2,2);
    plot(out,'o','Color',[0.6,0.6,0.8])
    x = quantile(k2',[0.25,0.75])';
    x = abs(x-[out,out]);
    hold on
    errorbar(1:20,out,x(1:20,1),x(1:20,2),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0,21])
    ylim([0 10e-11])
    set(gca,'ytick',0:2e-11:10e-11)
    box off
    title('PHS off')
    xlabel('Soil Layer')
    ylabel('Soil-to-root conductance (1/s)')
    
    
end






if ff(4)>0
    
    k2 = nan*smp(41:60,:);
    for i=41:60
        
        ix          = fsds>1;
        k2(i-40,ix) = qrootsink(i,ix)./min(189000,(smp(i,ix)+255000));
        
        ix1       = fsds>1&smp(i,:)<=-255000;
        k2(i-40,ix1) = 0;
        
        
    end
    
    xdk = figure;
    ix = fsds>1;
    subplot('position',[0.54, 0.12, 0.43, 0.8])
    h = boxplot(ksr(1:20,ix)');
    set(h(1:4,:),'Visible','off')
    set(h(7,:),'Visible','off')
    box off
    
    title('PHS on')
    xlabel('Soil Layer')
    text(18,6.5e-9,'(b)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:5:20)
    set(gca,'xticklabel',0:5:20)
    xlim([0 21])
    
    subplot('position',[0.08, 0.12, 0.43, 0.8])
    h = boxplot(k2');
    set(h(1:4,:),'Visible','off')
    set(h(7,:),'Visible','off')
    ylim([0, 7e-11])
    
    xlim([0,21])
    text(1.5,6.5e-11,'(a)','FontSize',14,'FontWeight','bold')
    
    text(18,6.5e-9,'(b)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:5:20)
    set(gca,'xticklabel',0:5:20)
    xlim([0 21])
    
    box off
    title('PHS off')
    xlabel('Soil Layer')
    ylabel('Soil-to-root conductance (1/s)')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(4)>1
        print(xdk,'../figs/fig3','-dpdf')
    end
    
    
end








if ff(5) >100
    xdk = figure;
    
    subplot('position',[0.07, 0.56, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    set(gca,'xticklabel',[])
    ylabel('Soil-root conductance')
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    set(gca,'xticklabel',[])
    
    subplot('position',[0.07, 0.12, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    xlabel('Soil Layer')
    ylabel('Soil-root conductance')

    subplot('position',[0.54, 0.12, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    xlabel('Soil Layer')
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(5)>1
        print(xdk,'../figs/fig5','-dpdf')
    end
end

