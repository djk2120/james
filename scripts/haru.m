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

ff = [0,0,0,0,2];


if ff(5) >0
    xdk = figure;
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = fsds>1;
    kon(1:20,ix) = ksr(1:20,ix);
    for ss = 21:40
        ix  = fsds>1&year>2001&abs(smp(ss,:)-vegwp(8,:))>1000;
        kon(ss,ix) = qrootsink(ss,ix)./(smp(ss,ix)-vegwp(8,ix));
    end
    
    for ss = 41:80    
        if ss>60
            ix          = year>2001&fsds>1;
        else
            ix          = fsds>1;
        end
        
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
    ylim([0 8e-11])
    text(17,6/7*8e-11,'(c)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    hold on
    ll = 61:80;
    plot(1:20,out(ll,1),'x')
    errorbar(1:20,out(ll,1),out(ll,2),out(ll,3),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlabel('Soil Layer')
    ylim([0 8e-11])
    text(17,6/7*8e-11,'(d)','FontSize',14,'FontWeight','bold')
    

    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(5)>1
        print(xdk,'../figs/fig3','-dpdf')
    end
    
    
    xdk = figure;
    bar(out(21:40,1))
    xlim([0 20.5])
    xlabel('Soil Layer')
    ylabel('Conductance (s-1)')
    title('PHS-on, TFE')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    if ff(5)>1
        print(xdk,'../figs/fig3a','-dpdf')
    end
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


if ff(2)>0
    %alternate conductance plot with bars instead of boxes
    
    k2 = nan*smp(41:60,:);
    for i=41:60
        
        ix          = fsds>1;
        k2(i-40,ix) = qrootsink(i,ix)./min(189000,(smp(i,ix)+255000));
        
        ix1       = fsds>1&smp(i,:)<=-255000;
        k2(i-40,ix1) = 0;
        
   
    end
    
    xdk = figure;
    ix = fsds>1;
    subplot('position',[0.54, 0.12, 0.43, 0.83])
    
    out  = median(ksr(1:20,ix),2);
    hold on
    bar(out,'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    x = quantile(ksr(1:20,ix)',[0.25,0.75])';
    x = abs(x-[out,out]);
    errorbar(1:20,out,x(1:20,1),x(1:20,2),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0,21])
    title('PHS on')
    xlabel('Soil Layer')
    text(18,6.5e-9,'(b)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.08, 0.12, 0.43, 0.83])
    bar(nanmedian(k2,2))
    out  = nanmedian(k2,2);
    hold on
    bar(out,'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    x = quantile(k2',[0.25,0.75])';
    x = abs(x-[out,out]);
    errorbar(1:20,out,x(1:20,1),x(1:20,2),'x','Color',[0.7 0.1 0.1],'Marker','none')
    xlim([0,21])
    text(1.5,6.5e-11,'(a)','FontSize',14,'FontWeight','bold')
    
    box off
    title('PHS off')
    xlabel('Soil Layer')
    ylabel('Soil-to-root conductance (1/s)')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(2)>1
        print(xdk,'../figs/fig3','-dpdf')
    end
    
    
end



if ff(1)>0
    %figure 2, water potential from a dry day
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



if ff(5) >100
    xdk = figure;
    
    subplot('position',[0.07, 0.12, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    xlabel('Soil Layer')
    ylabel('Soil-root conductance')
    
    subplot('position',[0.07, 0.56, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    set(gca,'xticklabel',[])
    ylabel('Soil-root conductance')
    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    xlabel('Soil Layer')
    
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    plot(1:20,(1:20)*10^-9)
    set(gca,'xticklabel',[])
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if ff(5)>1
        print(xdk,'../figs/fig5','-dpdf')
    end
end

