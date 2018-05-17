close all

dir  = '../data/mar6/';
dir2 = '../data/apr17/';

files = {...
    %    [dir,'BR-CAX_I1PTCLM50_r270.clm2.h1.2001-01-01-00000.nc'];...
    %    [dir,'BR-CAX_I1PTCLM50_r270_on_tfe.clm2.h1.2001-01-01-00000.nc'];...
    %    [dir,'BR-CAX_I1PTCLM50_r270_off_amb.clm2.h1.2001-01-01-00000.nc'];...
    %    [dir,'BR-CAX_I1PTCLM50_r270_off_tfe.clm2.h1.2001-01-01-00000.nc']};
    [dir2,'BR-CAX_I1PTCLM50_r270v3_phs_amb.clm2.h1.2001-01-01-00000.nc'];...
    [dir2,'BR-CAX_I1PTCLM50_r270v3_phs_tfe_60.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r270v3_sms_amb.clm2.h1.2001-01-01-00000.nc'];...
    [dir,'BR-CAX_I1PTCLM50_r270v3_sms_tfe_60.clm2.h1.2001-01-01-00000.nc']...
    };


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
    vard(3,:)  = [1,1,1,1];
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

dz = zs(2:end)-zs(1:end-1);

%************************************************************************
%------------------------------------------------------------------------

ff = [0,0,0,0,0,...
    0,0,0,0,0,...
    0,0,1,0,0,...
    0];


if ff(7)>0
    ix = year==2003&month>8&month<12;
    %create cumulative profile
    out = zeros(860,4);
    for ee=1:4
        x=180*sum(qrootsink((1:20)+(ee-1)*20,ix),2);
        dzi = round(100*dz);
        ss = 1;
        
        for i=1:860
            if i/100>zs(ss+1)
                ss = ss+1;
            end
            out(i,ee)=x(ss)/dzi(ss);
        end
    end
    
    %create cumulative timeseries
    ct = 0;
    out2 = zeros(8,5000);
    for w=[0,1]
        if ~w
            ix = year==2003&month>8&month<12;
        else
            ix = year==2003&month>1&month<5;
        end
        n  = sum(ix);
        for d = [0,1]
            if ~d
                sv = 1:4;
                %sv = 1:5;
            else
                sv=5:20;
                %sv = 6:20;
            end
            for i=1:4
                ct = ct+1;
                out2(ct,1:n) = 180*cumsum(sum(qrootsink(sv+(i-1)*20,ix)));
            end
        end
    end
    
    %first subplot
    xdk = figure;
    subplot('Position',[0.06,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(d)'};
    th = [150,400,300];
    tlab   = {'PHS-amb','PHS-tfe','SMS-amb','SMS-tfe'};
    for i=1:4
        plot(cumsum(flipud(out(:,i))),-8.6:0.01:-0.01,...
            'LineStyle',s{i},'Color',c(i,:),'LineWidth',2)
        hold on
    end
    set(gca,'yticklabel',9:-1:0)
    ylabel('Depth (m)')
    xlabel({'Cumulative Root Water Uptake (cm)'})
    legend(tlab,'location','Southeast')
    text(1,-0.5,'(a)','FontWeight','bold','FontSize',14)
    box off
    
    %other 3
    ix = year==2003&month>8&month<12;
    n  = sum(ix);
    b = 0.68:-0.29:0.1;
    ct = 0;
    
    for i=1:3
        subplot('Position',[0.56,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end
            set(gca,'xticklabel',[])
        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
            xlabel('Day of 2003')
        end
        %text(242,0.09*th(i),t{i},'FontWeight','bold')
        text(242,0.09*th(i),t{i})
        text(330,0.01*th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.51,0.39,'Cumulative Water (cm)',...
        'FontSize',11,'Rotation',90);
    text(0.5,0.97,'2003 Dry Season: Sept-Oct-Nov',...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if ff(7)>1
        print(xdk,'../figs3/fig7','-dpdf')
    end
    
end


if ff(8)>0
    ix = year==2003&month>1&month<5;
    
    %create cumulative profile
    out = zeros(860,4);
    for ee=1:4
        x=180*sum(qrootsink((1:20)+(ee-1)*20,ix),2);
        dzi = round(100*dz);
        ss = 1;
        
        for i=1:860
            if i/100>zs(ss+1)
                ss = ss+1;
            end
            out(i,ee)=x(ss)/dzi(ss);
        end
    end
    
    %create cumulative timeseries
    ct = 0;
    out2 = zeros(8,5000);
    n  = sum(ix);
    for d = [0,1]
        if ~d
            sv = 1:4;
            %sv = 1:5;
        else
            sv=5:20;
            %sv = 6:20;
        end
        for i=1:4
            ct = ct+1;
            out2(ct,1:n) = 180*cumsum(sum(qrootsink(sv+(i-1)*20,ix)));
        end
    end
    
    xdk = figure;
    subplot('Position',[0.06,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(d)'};
    th = [0.9*40,30-0.1*40,0.9*110];
    tlab   = {'PHS-amb','PHS-tfe','SMS-amb','SMS-tfe'};
    for i=1:4
        plot(cumsum(flipud(out(:,i))),-8.6:0.01:-0.01,...
            'LineStyle',s{i},'Color',c(i,:),'LineWidth',2)
        hold on
    end
    xlim([-12,26])
    set(gca,'yticklabel',9:-1:0)
    ylabel('Depth (m)')
    xlabel({'Cumulative Root Water Uptake (cm)'})
    legend(tlab,'location','Southeast')
    text(20.5,-0.5,'(a)','FontWeight','bold','FontSize',14)
    box off
    ix = year==2003&month>1&month<5;
    n  = sum(ix);
    b = 0.68:-0.29:0.1;
    ct = 0;
    for i=1:3
        subplot('Position',[0.56,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end
            set(gca,'xticklabel',[])
        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
            xlabel('Day of 2003')
        end
        xlim([31,121])
        text(34,th(i),t{i})
        text(90,th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    ylim([0,110])
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.51,0.39,'Cumulative Water (cm)',...
        'FontSize',11,'Rotation',90);
    text(0.5,0.97,'2003 Wet Season: Feb-Mar-Apr',...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if ff(8)>1
        print(xdk,'../figs3/fig8','-dpdf')
    end
    
    
end


if ff(6)>0
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>0;
    kon(1:40,:) = ksr(1:40,:);
    
    
    for ss = 41:80
        
        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-254900;
        kon(ss,ix1) = 0;
    end
    
    xdk = figure;
    mm = 3;
    ix = year==2003;
    g  = findgroups(mcsec(ix));
    subplot(1,2,1)
    tv = 0.25:0.5:24;
    x=splitapply(@nanmean,kon(43,ix),g);
    plot(tv,x,'k','LineWidth',2)
    hold on
    x=splitapply(@nanmean,kon(63,ix),g);
    plot(tv,x,'k:','LineWidth',2)
    set(gca,'xtick',6:6:18)
    xlim([6,18])
    title('SMS')
    xlabel('Hour of Day')
    ylabel({'Average (implied)';'Hydraulic Conductance (s^{-1})'})
    text(16.5,4.7e-11,'(a)','FontWeight','bold','FontSize',14)
    
    subplot(1,2,2)
    x=splitapply(@nanmean,kon(3,ix),g);
    plot(tv,x,'k','LineWidth',2)
    hold on
    x=splitapply(@nanmean,kon(23,ix),g);
    plot(tv,x,'k:','LineWidth',2)
    set(gca,'xtick',6:6:18)
    xlim([6,18])
    ylim([0,5e-9])
    title('PHS')
    xlabel('Hour of Day')
    ylabel({'Average (modeled)';'Hydraulic Conductance (s^{-1})'})

    text(16.5,4.7e-9,'(b)','FontWeight','bold','FontSize',14)

     
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,3];
    xdk.PaperSize = [7,3];
    xdk.PaperPosition = [0,0,7,3];
    
    if ff(6)>1
        print(xdk,'../figs3/fig6','-dpdf')
    end
    
    
end


if ff(5)>0
    
    out = zeros(80,3);
    kon = nan*smp(1:80,:);
    ix  = year==2003&fctr(1,:)>0;
    kon(1:40,:) = ksr(1:40,:);
    
    
    for ss = 41:80
        
        kon(ss,ix)  = qrootsink(ss,ix)./min(189000,(smp(ss,ix)+255000));
        ix1         = ix&smp(ss,:)<=-254900;
        kon(ss,ix1) = 0;
    end
 
    ix = year==2003&(month==7|month==8|month==9);
    g  = findgroups(doy(ix));
    
    xdk = figure;
    
    dv = unique(doy(ix));
    subplot(3,1,2)
    plot(dv,splitapply(@mean,kon(3,ix),g),'k','LineWidth',2)
    ylabel({'Daily Mean (Modeled)';'Conductance (1/s)'})
    xlim([180,274])
    text(182,1e-9,'PHSamb')
    
    subplot(3,1,1)
    plot(dv,splitapply(@nanmean,kon(43,ix),g),'k','LineWidth',2)
    title('Jul-Aug-Sept 2003','FontSize',14)
    ylim([0,2e-11])
    xlim([180,274])
    text(182,(1/3)*10^-11,'SMSamb')
    ylabel({'Daily Mean (Implied)';'Conductance (1/s)'})
        set(gca,'ytick',0:1e-11:2e-11)
    
    
    subplot(3,1,3)
    bar(dv,1800*splitapply(@sum,prec(ix),g))
    xlim([180,274])
    xlabel('Day of Year')
    ylabel('Rain (mm)')

    
     xdk.Units = 'inches';
    xdk.Position = [2,2,7,6];
    xdk.PaperSize = [7,6];
    xdk.PaperPosition = [0,0,7,6];
    
    if ff(5)>1
        print(xdk,'../figs3/suppcond','-dpdf')
    end
    
end

if ff(13)>0
    ix1 = year==2003&month>1&month<5;
    ixd = mcsec>=diurn(tt)&mcsec<=diurn(tt+3);
    
    
    subplot(2,1,1)
    x1 = mean(qrootsink(1:20,ix1&ixd),2);
    
    x2 = mean(qrootsink(21:40,ix1&ixd),2);
    
    bar([x1,x2])
    
    figure
    subplot(1,2,1)
    bar(mean(smp(1:20,ix1)/101972,2))
    ylim([-0.3,0])
    subplot(1,2,2)
    bar(mean(smp(21:40,ix1)/101972,2))
    ylim([-0.3,0])
    
end

if ff(14)>0
    xdk = figure;
    c='PHSambPHStfeSMSambSMStfe';
    x= [-0.19,0.95*-2.5,-0.19,0.95*-2.5];
    tt=25;
    ixd = mcsec>=diurn(tt)&mcsec<=diurn(tt+3);
    ix = year==2003&ixd;
    s  = [1,3,2,4];
    for i = 1:4
        subplot(2,2,i)
        plot(smp(3+(i-1)*20,ix)/101972,qrootsink(3+(i-1)*20,ix),'.')
        if i==1||i==3
            xlim([-0.2,0])
            %ylim([-1e-5,5e-5])
            ylim([-5e-6,10e-5])
            %xlim([-2.5,0])
            ylabel({'Root Water Uptake';'(mm/s)'})
        else
            xlim([-2.5,0])
            ylim([-5e-6,10e-5])
        end
        if i>2
            xlabel('Soil Potential')
        end
        title(c((1:6)+(i-1)*6))
        
        grid on
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,6];
    xdk.PaperSize = [7,6];
    xdk.PaperPosition = [0,0,7,6];
    
    if ff(14)>1
        print(xdk,'../figs3/supprwu','-dpdf')
    end
    
end


if ff(15)>0
   
    tstr = {'PHSamb','PHStfe','SMSamb','SMStfe'}
    xdk = figure;
    
    for i=1:4
        
        ymin = [0,-0.1,0,0];
        ymax = [0.05,0.5,3,3];
        
        if i==1||i==2
        x = smp(3+(i-1)*20,:)-vegwp(4+(i-1)*4,:);
        else
            x = smp(3+(i-1)*20,:)+255000;
        end
        
        tv = 0.25:0.5:24;
        ix = year==2003; 
        g = findgroups(mcsec(ix));
        subplot(2,2,i)
        plot(tv,splitapply(@mean,x(ix),g)/101972,'k','LineWidth',1.5)
        xlim([6,18])
        set(gca,'xtick',6:3:18)
        ylim([ymin(i),ymax(i)])
          title(tstr{i})  
        if i>2
            xlabel('Hour of Day')
        end
        if i==1||i==3
            ylabel({'Average  \Delta\psi (MPa)'})
        end
        
    end
    

        
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,6];
    xdk.PaperSize = [7,6];
    xdk.PaperPosition = [0,0,7,6];
    
    if ff(15)>1
        print(xdk,'../figs3/supppsi','-dpdf')
    end
    
    
    
end

