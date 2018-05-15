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
    0,0,1];


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
    subplot('Position',[0.06,0.13,0.42,0.83])
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
    b = 0.71:-0.29:0.1;
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
        text(242,0.09*th(i),t{i},'FontWeight','bold')
        text(330,0.01*th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.51,0.39,'Cumulative Water (cm)',...
        'FontSize',11,'Rotation',90);

    
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
    subplot('Position',[0.06,0.13,0.42,0.83])
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
    b = 0.71:-0.29:0.1;
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
        text(34,th(i),t{i},'FontWeight','bold')
        text(90,th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    ylim([0,110])
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.51,0.39,'Cumulative Water (cm)',...
        'FontSize',11,'Rotation',90);
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if ff(7)>1
        print(xdk,'../figs3/fig7','-dpdf')
    end
    
    
end


