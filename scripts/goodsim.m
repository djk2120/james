close all

dir  = '../data/mar6/';
dir2 = '../data/apr17/';
dir3 = '/Users/kennedy/Desktop/james/goodsim/may18/';

files = {...
    [dir3,'BR-CAX_I1PTCLM50_r270v5_phs_amb_bf100_n28.clm2.h1.2001-01-01-00000.nc'];...
    [dir3,'BR-CAX_I1PTCLM50_r270v5_phs_tfe60_bf100_n28.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_amb_bf100_n28_kr1e-10.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_tfe60_bf100_n28_kr1e-10.clm2.h1.2001-01-01-00000.nc']...
    [dir3,'BR-CAX_I1PTCLM50_r270v8_phs_amb_bf100_n28_hk2.clm2.h1.2001-01-01-00000.nc'];...
    [dir3,'BR-CAX_I1PTCLM50_r270v8_phs_tfe60_bf100_n28_hk2.clm2.h1.2001-01-01-00000.nc']...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_amb_bf100_n28_kr1e-9.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_tfe60_bf100_n28_kr1e-9.clm2.h1.2001-01-01-00000.nc']...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_amb_bf100_n28_kr3.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v6_phs_tfe60_bf100_n28_kr3.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v4_phs_amb_bf100.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v4_phs_tfe_60_bf100.clm2.h1.2001-01-01-00000.nc'];...
    %[dir3,'BR-CAX_I1PTCLM50_r270v3_phs_tfe_60.clm2.h1.2001-01-01-00000.nc'];...
    %[dir,'BR-CAX_I1PTCLM50_r270v3_sms_amb.clm2.h1.2001-01-01-00000.nc'];...
    %[dir,'BR-CAX_I1PTCLM50_r270v3_sms_tfe_60.clm2.h1.2001-01-01-00000.nc']...
    };


a=ncinfo(files{1});

offset = 10;
ns     = 20;
nx     = length(files);

if ~exist('fctr','var')
    a          = getvars( files{1} , offset, ns, 1);
    varlist    = {'FCTR','FPSN','BTRAN','VEGWP','SMP','QROOTSINK',...
        ...              1      2      3       4        5     6
        'FCEV','FSH','KSR','FGEV','H2OSOI','ELAI','QOVER','QDRAI'};
    %     7      8      9     10      11     12     13      14
    vard       = ones(length(varlist),length(files));
    vard(3,:)  = [1,1,1,1];
    vard(4,:)  = [4,4,4,4];
    vard(5,:)  = ns;
    vard(6,:)  = ns;
    vard(9,:)  = [ns,ns,ns,ns];
    vard(11,:) = ns;
    vard(12,:) = [1,0,0,0];
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
    0,0,0,0,1,...
    0,0,0,0,0,...
    0];

if ff(10)>0
    
    params = cell(7,4);
    params(1,:) = {'rootprof_beta',0.95,0.98,0.993};
    params(2,:) = {'ck',2.95,3.95,5.45};
    params(3,:) = {'kmax',2e-8,4e-8,8e-8};
    params(4,:) = {'krmax',2e-9,6e-9,1.8e-8};
    params(5,:) = {'p50_34',-178451,-229437,-280423};
    params(6,:) = {'medlynslope',6,7,nan};
    params(7,:) = {'p50_+',50986,0,nan};
    
    n = 0;
    og = '../data/k4g7.nc';
    for a = 1:3
    for b = 1:3
    for c = 1:3
    for d = 1:3
    for e = 1:3
    for f = 1:2
    for g = 1:2
        n = n+1;
        p = [a,b,c,d,e,f,g];
        
        if n<10
            nstr = ['00',num2str(n)];
        elseif n<100
            nstr = ['0',num2str(n)];
        else
            nstr = num2str(n);
        end
        
        if n<2
            filename = ['n',nstr,'.nc'];
            fullpath = ['pfiles/',filename];
            cmd = ['cp ',og,' ',fullpath];
            system(cmd);
            
            for i = [1:4,6]
                x = ncread(fullpath,params{i,1});
                x(:) = params{i,1+p(i)};
                ncwrite(fullpath,params{i,1},x)
            end
            p50 = repmat(params{5,1+p(5)},79,4);
            p50(:,1:2) = p50(:,1:2)+params{7,1+p(7)};
            disp(p50)
            
        end
        
        
    end
    end
    end
    end
    end
    end
    end
    
        
    
end

if ff(11)>0
    ix = year>2001;
    g = month+(year-2002)*12;
    
    out = 1800*4e-7*splitapply(@sum,fctr(:,ix)',g(ix)')./repmat(eomday(2001,1:12)',2,4);
    %out = splitapply(@mean,fctr(:,ix)',g(ix)');
    
    pp = nan(24,2);
    mx = cumsum(repmat(eomday(2001,1:12),1,2));
    ll = 0;
    oneto = (1:730)';
    ix1 = p(:,1)>0;
    ix2 = p(:,2)>0;
    for i=1:24
        ix = oneto>ll&oneto<=mx(i);
        pp(i,1) = mean(p(ix&ix1,1));
        pp(i,2) = mean(p(ix&ix2,2));
        ll = mx(i);
    end
    
    
    figure
    plot(pp(1:12,2)./pp(1:12,1),'kx')
    hold on
    plot(pp(13:24,2)./pp(13:24,1),'rx')
        set(gca,'xtick',0:3:12)
            grid on
    ylim([0,1.3])
            
    figure
    out = splitapply(@mean,fctr',(year'-2001)*12+month');
    plot(out(13:24,2)./out(13:24,1),'kx')
    hold on
    plot(out(25:36,2)./out(25:36,1),'rx')
    grid on
    ylim([0,1.3])    
    
        figure
    out = splitapply(@mean,fctr',(year'-2001)*12+month');
    plot(out(13:24,4)./out(13:24,2),'kx')
    hold on
    plot(out(25:36,4)./out(25:36,2),'rx')
    grid on
    ylim([0,1.3])    
    
    figure
    out = splitapply(@mean,vpd,(year-2001)*12+month);
    plot(reshape(out,12,3))
    xlim([0,13])
    grid on
    
    figure
    ix = mcsec==diurn(10);
    g = (year-2001)*12+month;
    out = splitapply(@mean,vegwp(16,ix),g(ix));
    plot(reshape(out,12,3))
    xlim([0,13])
    grid on
    
end

if ff(15)>0
    ct = 0;
    for yy=2002:2003
        for mm=1:12
            ct = ct+1;
            subplot(4,6,ct)
            ix = year==yy&month==mm;
            plot(mean(smp(21:40,ix)/101972,2))
            hold on
            plot(mean(smp(61:80,ix)/101972,2))
            xlim([0,21])
            ylim([-1,0])
            title(mm)
        end
    end
    
    figure
    plot(vegwp(8,year>2001&mcsec==diurn(10)))
    hold on
    plot(vegwp(16,year>2001&mcsec==diurn(10)))
    
    

    
    
end

if ff(1)>0
    
    ix = year==2003;
    tv = doy(ix);
    subplot(1,2,1)
    
    hold on
    plot(tv,180*cumsum(sum(qrootsink(61:67,ix))))
    plot(tv,180*cumsum(sum(qrootsink(21:27,ix))))
    xlim([0,365])
    set(gca,'xtick',0:90:365)
title('above 68cm')
    ylabel('Cumulative Water Extraction (cm)')
    xlabel('Day of 2003')
     subplot(1,2,2)
    
    hold on
    plot(tv,180*cumsum(sum(qrootsink(68:80,ix))))
    plot(tv,180*cumsum(sum(qrootsink(28:40,ix))))
    xlim([0,365])
        set(gca,'xtick',0:90:365)
        ylabel('Cumulative Water Extraction (cm)')
        xlabel('Day of 2003')
        title('below 68cm')
    
end


if ff(2) >0
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    a = csvread('../goodsim/control_sm.csv');
    a(a==0) = nan;
    b = csvread('../goodsim/tfe_sm.csv');
    b(b==0) = nan;
    
    p = {'PHSamb \theta_{s}=0.28','PHStfe \theta_{s}=0.28',...
        'PHSamb \theta_{s}=0.42','PHStfe \theta_{s}=0.42'};
    xdk = figure;
    
    for i=1:4
        subplot(2,2,i)
        plot(tv,100*h2osoi(17+(i-1)*20,:))
        xlim([0,1095])
        set(gca,'xtick',0:365/2:1095)
        set(gca,'xticklabel',{0,'',365,'',730,'',1095})
        ylim([5,35])
        xlabel('Day')
        ylabel('Volumetric Soil Water')
        title(p{i})
        hold on
        if i==1||i==3
            plot(a(:,1),a(:,8),'rx')
        else
            plot(b(:,1),b(:,8),'rx')
        end
        if i==2
            legend({'Model','Obs'},'Location','Northeast')
        end
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5,0.97,['Depth = ','5m'],...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    fout = ['../goodsim/figs/soilwater_','5m','_bf100_n28vs42'];
    mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
    
    if ff(2)>1
        print(xdk,fout,'-dpdf')
        system(mkjpg)
    end
    
    
    
    
end

if ff(3)>0
    
    mstr = 'janfebmaraprmayjunjulaugsepoctnovdec'    ;
    
    a = csvread('../goodsim/control_sm.csv');
    a(a==0) = nan;
    b = csvread('../goodsim/tfe_sm.csv');
    b(b==0) = nan;
    
    


    ix = year>2000;
    g = (year-2001)*365+doy;
    out = 100*splitapply(@mean,h2osoi(:,ix)',g(ix)')';
    
    ll = [7,9,12,14,16,17];
    
    figure
    for i=3:8
        subplot(2,3,i-2)
        plot(a(:,1),a(:,i),'rx')
        hold on
        plot(out(ll(i-2),:))
        plot(out(ll(i-2)+40,:))
        ylim([5,35])
    end
    
    figure
    for i=3:8
        subplot(2,3,i-2)
        plot(b(:,1),b(:,i),'rx')
        hold on
        plot(out(ll(i-2)+20,:))
        plot(out(ll(i-2)+60,:))
        ylim([5,35])
    end
    
end

if ff(4)>0
    i = 0;
    for yy=2001:2003
        for mm=1:2:12
            i = i+1;
            subplot(3,6,i)
            ix = year==yy&month==mm;
            plot(mean(smp(21:40,ix),2)/101972)
            hold on
            plot(mean(smp(61:80,ix),2)/101972)
            ylim([-1,0])
        end
    end
    
    figure
    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,(fctr)',g');
    
    plot(cumsum(out(:,2)))
    hold on
    plot(cumsum(out(:,4)))
    
    figure
    plot(vegwp(8,mcsec==diurn(10)))
    hold on
    plot(vegwp(16,mcsec==diurn(10)))
    
    figure
    
    out1=(splitapply(@mean,ksr(25,:),g));
    out2=(splitapply(@mean,ksr(65,:),g));
    hold on
    plot(out2./out1,'.')
    set(gca,'xtick',0:365/2:1095)
    ylim([0,4])
    
    
end

if ff(5)>0
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    
    xdk = figure;
    for i=1:5
        subplot(2,5,i)
        
        plot(tv,1800*cumsum(qrootsink(i+35,:)))
        
        xlim([0,1095])
        set(gca,'xtick',0:365:1095)
        
        
        if i==1
            ylabel({'Cumulative Water','Extracted (mm)'})
            title(['Soil layer ',num2str(i+15)])
        else
            title(num2str(i+15))
        end
        
        subplot(2,5,i+5)
        plot(tv,1800*cumsum(qrootsink(i+35,:)))
        xlim([0,1095])
        set(gca,'xtick',0:365/2:1095)
        set(gca,'xticklabel',{0,'',365,'',730,'',1095})
        xlabel('Day')
        
        if i==1
            ylabel({'Cumulative Water','Extracted (mm)'})
        end
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.85,0.85,{'PHSamb, bf*100';'saturation = 28%'},'FontWeight','bold')
    text(0.85,0.15,{'PHSamb, bf*100';'saturation = 42%'},'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,9,4];
    xdk.PaperSize = [9,4];
    xdk.PaperPosition = [0,0,9,4];
    
    
    fout = ['../goodsim/figs/extraction42v28'];
    mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
    
    if ff(5)>1
        print(xdk,fout,'-dpdf')
        system(mkjpg)
    end
    
    
end

if ff(6)>0
    ix = year>2001;
    g  = (year-2002)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr(:,ix)',g(ix)');
    
    oneto = (1:730)';

    
    ss = [1,3,2,4];
    for i=1:4
        subplot(2,2,ss(i))
        
        if i==1||i==3
            ix = p(:,1)>0&oneto>365;
            plot(oneto(ix),out(ix,i),'.')
            hold on
            plot(oneto(ix),p(ix,1),'.')
        else
            ix = p(:,2)>0&oneto>365;
            plot(oneto(ix),out(ix,i),'.')
            hold on
            plot(oneto(ix),p(ix,2),'.')
        end
    end
    
    if 1==2
    figure
    subplot(2,2,1)
    plot([0,6],[0,6],'k')
    hold on
    plot(out(ix,1),p(ix,1),'.')
    xlabel('Model')
    ylabel('OBS')
    xlim([0,6])
    ylim([0,6])
    subplot(2,2,2)
    plot([0,6],[0,6],'k')
    hold on
    plot(out(ix,2),p(ix,2),'.')
    xlabel('Model')
    ylabel('OBS')
    xlim([0,6])
    ylim([0,6])
        subplot(2,2,3)
    plot([0,6],[0,6],'k')
    hold on
    plot(out(ix,3),p(ix,1),'.')
    xlabel('Model')
    ylabel('OBS')
    xlim([0,6])
    ylim([0,6])
    subplot(2,2,4)
    plot([0,6],[0,6],'k')
    hold on
    plot(out(ix,4),p(ix,2),'.')
    xlabel('Model')
    ylabel('OBS')
    xlim([0,6])
    ylim([0,6])
    end
end

if ff(7)>0
    
    ix = year==2003&mcsec==diurn(26)&month>1&month<5;
    subplot(2,2,1)
    hold on
    plot(vegwp(1,ix)/101972)
    plot(vegwp(9,ix)/101972)
    ylim([-3,0])
    subplot(2,2,2)
    hold on
    plot(vegwp(5,ix)/101972)
    plot(vegwp(13,ix)/101972)
    ylim([-3,0])
    
    
    ix = year==2003&month>1&month<5;
    g = findgroups(doy(ix));
    out = splitapply(@mean,fctr(:,ix)',g');
    
    subplot(2,2,3)
    plot(out(:,1))
    hold on
plot(out(:,3))
    subplot(2,2,4)
    plot(out(:,2))
    hold on
plot(out(:,4))

    
    
end

if ff(8)>0
    t = 0.05:0.01:0.42;
    tsat = 0.42;
    psat = -47/101.972;

    b=6;
    ksat = 0.03*3600;
    p = psat*(t/tsat).^-b;
    c = 2*b+3;
    k = ksat*(t/tsat).^c;
    
    plot(log10(-p),log10(k))
    set(gca,'xtick',-1:1)
    xlim([-1.5,1.5])
    ylim([-2.5,3.5])
    hold on
    grid on
    
    figure
    plot(log10(-p),t)
    ylim([0,0.4])
    xlim([0.5,3.5])
    grid on
    
end

if ff(9)>0
    ix = year>2001;
    g = (year(ix)-2002)*12+month(ix);
    out = 4e-7*86400*splitapply(@mean,fctr(:,ix)',g');
    xv = (cumsum([0,eomday(2001,1:11)])+cumsum(eomday(2001,1:12)))/2;
    xv = ([xv,365+xv])/365+2002;
    
    xdk = figure;
    for i=1:4
        subplot(2,2,i)

    set(gca,'xtick',2002:0.5:2004)
    oneto = (1:730)/365+2002;

    hold on
    if i==1||i==3
                ix = p(:,1)>0;
    plot(oneto(ix),p(ix,1),'rx')
    else
        ix = p(:,2)>0;
    plot(oneto(ix),p(ix,2),'rx')
    end
        plot(xv,out(:,i),'LineWidth',2)
    ylabel('Transpiration (mm/d)')
    end
    

    
    figure
    ix = year>2001;
    g  = doy(ix)+(year(ix)-2002)*365;
    out = 1800*4e-7*splitapply(@sum,fctr(:,ix)',g');
    
    oneto = (1:730)/365+2002;
    for i=1:4
        subplot(2,2,i)
        if i==1||i==3
                        ix = p(:,1)>0;
                plot(oneto(ix),out(ix,i)-p(ix,1),'.')
        else
            ix = p(:,2)>0;
                plot(oneto(ix),out(ix,i)-p(ix,2),'.')
        end

    ylabel('Transpiration Model-OBS (mm/d)')
    end
    
end

if ff(12)>0
    
    somestrings = {'50cm','1m','2m','3m','4m','5m'};
    
    
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    
    a = csvread('../goodsim/control_sm.csv');
    a(a==0) = nan;
    b = csvread('../goodsim/tfe_sm.csv');
    b(b==0) = nan;
    
    %p = {'PHSamb','PHStfe','SMSamb','SMStfe'};
    p = {'PHSamb','PHStfe','PHSamb bf*100','PHStfe bf*100'};
    
    oneto = 1:21;
    
    tt = zeros(6,1);
    i = 0;
    for ll = [0.5,1:5]
        i = i+1;
        tt(i) = min(oneto(zs>ll))-1;
    end
    
    for dd=1:6
        if 1==1
            xdk = figure;
            
            
            
            for i=1:4
                subplot(2,2,i)
                
                if i==1||i==3
                    targ = a;
                else
                    targ = b;
                end
                
                
                plot(tv,100*h2osoi(tt(dd)+(i-1)*20,:))
                hold on
                plot(targ(:,1),targ(:,dd+1),'rx')
                ylim([5,45])
                xlim([0,1095])
                set(gca,'xtick',0:365/2:1095)
                set(gca,'xticklabel',{0,'',365,'',730,'',1095})
                xlabel('Day')
                ylabel('Volumetric Soil Water')
                title(p{i})
                
                if i==2
                    legend({'Model','Obs'},'Location','Northeast')
                end
            end
            
            ax1 = axes('Position',[0 0 1 1],'Visible','off');
            text(0.5,0.97,['Depth = ',somestrings{dd}],...
                'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
            
            
            xdk.Units = 'inches';
            xdk.Position = [2,2,7,5];
            xdk.PaperSize = [7,5];
            xdk.PaperPosition = [0,0,7,5];
            
            %fout = ['../goodsim/figs/soilwater_',somestrings{dd}];
            fout = ['../goodsim/figs/soilwater_',somestrings{dd},'_bf100'];
            mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
            
            if ff(12)>0
                print(xdk,fout,'-dpdf')
                system(mkjpg)
            end
        end
        
    end
end


if ff(13)>0
    
    xdk = figure;
    
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    
    subplot(2,2,1)
    a = ncread(files{1},'QOVER');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    hold on
    a = ncread(files{2},'QOVER');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    ylim([0,2400])
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    
    legend('AMB','TFE','location','Northwest')
    xlabel('Day')
    ylabel('Cumulative QOVER (mm)')
    title('PHS')
    box on
    
    subplot(2,2,2)
    hold on
    a = ncread(files{3},'QOVER');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    a = ncread(files{4},'QOVER');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    ylim([0,2400])
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    
    legend('AMB','TFE','location','Northwest')
    xlabel('Day')
    ylabel('Cumulative QOVER (mm)')
    title('PHS, bf*100')
    box on
    
    subplot(2,2,3)
    a = ncread(files{1},'QDRAI');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    hold on
    a = ncread(files{2},'QDRAI');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    ylim([0,2400])
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    
    legend('AMB','TFE','location','Northwest')
    xlabel('Day')
    ylabel('Cumulative QDRAI (mm)')
    title('PHS')
    box on
    
    subplot(2,2,4)
    hold on
    a = ncread(files{3},'QDRAI');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    a = ncread(files{4},'QDRAI');
    a = a(1+offset:end);
    plot(tv,cumsum(1800*a))
    ylim([0,2400])
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    legend('AMB','TFE','location','Northwest')
    xlabel('Day')
    ylabel('Cumulative QDRAI (mm)')
    title('PHS, bf*100')
    box on
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    fout = ['../goodsim/figs/drainage'];
    mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
    
    if ff(13)>0
        print(xdk,fout,'-dpdf')
        system(mkjpg)
    end
    
    
end




if ff(14)>0
    
    xdk = figure;
    tv = doy+mcsec/max(diurn)-1+(year-2001)*365;
    a = ncread(files{1},'QOVER');
    plot(tv,cumsum(1800*a(1+offset:end)),'-','LineWidth',2)
    a = ncread(files{2},'QOVER');
    hold on
    plot(tv,cumsum(1800*a(1+offset:end)),'-','LineWidth',2)
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    xlabel('Day')
    ylabel('Cumulative QOVER (mm)')
    title('Baseflow x100')
    legend('AMB','TFE','location','southeast')
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    xdk = figure;
    tv = doy+mcsec/max(diurn)-1+(year-2001)*365;
    a = ncread(files{1},'QDRAI');
    plot(tv,cumsum(1800*a(1+offset:end)),'-','LineWidth',2)
    a = ncread(files{2},'QDRAI');
    hold on
    plot(tv,cumsum(1800*a(1+offset:end)),'-','LineWidth',2)
    xlim([0,1095])
    set(gca,'xtick',0:365/2:1095)
    set(gca,'xticklabel',{0,'',365,'',730,'',1095})
    xlabel('Day')
    ylabel('Cumulative QDRAI (mm)')
    title('Baseflow x100')
    legend('AMB','TFE','location','southeast')
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    
    
    
    
end


