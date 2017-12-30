close all

dir = '../data/aug25/';
files = {...
    [dir,'BR-CAX_I1PTCLM50_r251_k5g7.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r251tfe_k5g7.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r251off_k5g7.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r251offtfe_k5g7.clm2.h1.2001-01-01-00000.nc']};

a=ncinfo(files{1});

offset = 10;
ns     = 20;
nx     = length(files);

if ~exist('fctr','var')
    a          = getvars( files{1} , offset, ns, 1);
    varlist    = {'FCTR','FPSN','BTRAN','VEGWP','SMP','QROOTSINK',...
        ...              1      2      3       4        5     6
        'FCEV','FSH','KSR'};
    %     7      8      9     10
    vard       = ones(length(varlist),length(files));
    vard(3,:)  = [0,0,1,1];
    vard(4,:)  = [4,4,0,0];
    vard(5,:)  = ns;
    vard(6,:)  = ns;
    vard(9,:) = [ns,0,0,0];
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

%************************************************************************
%------------------------------------------------------------------------

ff = [2];

if ff(1)>0
    dd=2;
    mm=11;
    yy=2002;
    
    xdk = figure;
    x = 0.25:0.5:24;
    xf = 1/101972;
    subplot('position',[0.08, 0.12, 0.43, 0.83])
    plot(x,xf*vegwp(1:4,day==dd&month==mm&year==yy)')
    xlim([0 24])
    ylim([-2.5 0])
    xlabel('Hour')
    ylabel('Water Potential (MPa)')
    title('AMB')
    text(1.5,-2.3,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:6:24)
    subplot('position',[0.54, 0.12, 0.43, 0.83])
    plot(x,xf*vegwp(5:8,day==dd&month==mm&year==yy)')
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
    if ff(1)>1
        print(xdk,'../figs/fig2','-dpdf')
    end
    
    
    % tfe pressure drop is about 2x larger soil-to-root as compared to
    % root-to-leaf
    if 1==2
        figure
        plot(vegwp(1,day==dd&month==mm&year==yy)-vegwp(4,day==dd&month==mm&year==yy)...
            -(vegwp(5,day==dd&month==mm&year==yy)-vegwp(8,day==dd&month==mm&year==yy)))
        hold on
        x1 = vegwp(4,day==dd&month==mm&year==yy);
        x2 = vegwp(8,day==dd&month==mm&year==yy);
        plot(x2-max(x2)-(x1-max(x1)))
    end
    
end
