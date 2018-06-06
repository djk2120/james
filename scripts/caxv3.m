close all

dir  = '../data/mar6/';
dir2 = '../data/apr17/';
dir3 = '../data/may31/';

%files = {...
%    [dir3,'BR-CAX_r270_amb_n374.nc'];...
%    [dir3,'BR-CAX_r270_tfe_n374.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r270v3_sms_amb.clm2.h1.2001-01-01-00000.nc'];...
%    [dir,'BR-CAX_I1PTCLM50_r270v3_sms_tfe_60.clm2.h1.2001-01-01-00000.nc']...
%    };
files = {...
    [dir3,'BR-CAX_I1PTCLM50_r270v9_phs_amb.clm2.h1.2001-01-01-00000.nc'];...
    [dir3,'BR-CAX_I1PTCLM50_r270v9_phs_tfe.clm2.h1.2001-01-01-00000.nc'];...
    [dir3,'BR-CAX_I1PTCLM50_r270v9_sms_amb.clm2.h1.2001-01-01-00000.nc'];...
    [dir3,'BR-CAX_I1PTCLM50_r270v9_sms_tfe.clm2.h1.2001-01-01-00000.nc']...
    };



a=ncinfo(files{1});

offset = 10;
ns     = 20;
nx     = length(files);

if ~exist('fctr','var')
    a          = getvars( files{3} , offset, ns, 1);
    varlist    = {'FCTR','FPSN','BTRAN','VEGWP','SMP','QROOTSINK',...
        ...          1      2      3       4        5     6
        'FCEV','FSH','KSR','ELAI','FGEV','H2OSOI'};
    %     7      8      9     10      11     12     13      14
    vard       = ones(length(varlist),length(files));
    vard(3,:)  = [1,1,1,1];
    vard(4,:)  = [4,4,0,0];
    vard(5,:)  = ns;
    vard(6,:)  = ns;
    vard(9,:)  = [ns,ns,0,0];
    vard(10,:) = [0,0,1,0];
    vard(12,:) = ns;
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

zr = 1-0.98^(zs(2)*100);
for i=2:20
    zr(i) = 1-sum(zr)-0.98^(zs(i+1)*100);
end


if ~exist('p','var')
f = '../data/daily_sapflow_control.txt';
a=csvread(f,1,0);
p = a(:,2);
d = a(:,1);
f = '../data/daily_sapflow_drought.txt';
a=csvread(f,1,0);
p = [p,a(:,2)];
p = [zeros(365,2);p];

dz = zs(2:end)-zs(1:end-1);
end
%************************************************************************
%------------------------------------------------------------------------

hh = zeros(200,1);
gg = zeros(200,1);
ff = zeros(200,1);

hh(1:5)  = [1,0,0,0,0];
hh(6:10) = [0,0,0,0,0];

if hh(10)>0
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    g  = (year(ix)-2001)*365+doy(ix);
    subplot(2,1,1)
    plot(splitapply(@mean,vpd(ix),g),'.')
    set(gca,'xtick',cumsum([0,eomday(2001,1:12),eomday(2001,1:12),eomday(2001,1:12)]))
    set(gca,'xticklabel',[0:12,1:12,1:12])
    xlim([0,1095])
    grid on
    
    
    subplot(2,1,2)
    plot(splitapply(@mean,fsds(ix),g),'.')
    set(gca,'xtick',cumsum([0,eomday(2001,1:12),eomday(2001,1:12),eomday(2001,1:12)]))
    set(gca,'xticklabel',[0:12,1:12,1:12])
    xlim([0,1095])
    grid on
    
    
    
    
end

if hh(9)>0
    xv = 0.5:730;
    xdk = figure;
    
    
    years =cell(1,25);
    for i=0:24
        if i==6
            years(i+1)={'2002'};
        elseif i==18
            years(i+1)={'2003'};
        end
    end
    
    mls =cell(1,25);
    mm = -1;
    for i=0:24
        if mod(i,2)~=0
            mm = mm+2;
            if mm>12
                mm=1;
            end
            mls(i+1)={num2str(mm)};
        end
    end
    
    
    xx = [0.09,0.54,0.09,0.54];
    yy = [0.78,0.78,0.55,0.55];
    w  = 0.44;
    h  = 0.17;
    cc = [0,0,0.7,0.7];
    ss = '(a)(b)(c)(d)';
    for i=1:4
        s=subplot('Position',[xx(i),yy(i),w,h]);

        ix = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28);
        g = (year(ix)-2002)*365+doy(ix);
        plot(xv,splitapply(@mean,btran(i,ix),g),'.','color',ones(1,3)*cc(i))
                hold on
                plot([365,365],[0,1],'Color',[0.3,0.3,0.3])

        xlim([0,730])
        ylim([0,1])
        set(gca,'xtick',cumsum([0,eomday(2001,1:12),eomday(2001,1:12)]))
        set(gca,'xticklabel',years)
        grid on
        if i==1
            title('AMB')
            l = legend('PHS');
            l.Position =[0.32    0.7944    0.11    0.04];
        elseif i==2
            title('TFE')
        elseif i==3
            ll=legend('SMS');
            ll.Position =[0.32    0.5644    0.11    0.04];
        end

        if i==2||i==4
            set(gca,'yticklabel',[])
        end
        text(10,0.2,ss((1:3)+(i-1)*3),'FontSize',14,'FontWeight','bold')
        box off
    end
    
    xf = 86400*1e-6*12; %umol/m2/s to g/m2/d
    yy = [0.32,0.32,0.09,0.09];
    ss = '(e)(f)(g)(h)';
    for i=1:4
        s=subplot('Position',[xx(i),yy(i),w,h]);
                plot([365,365],[0,10],'Color',[0.3,0.3,0.3])
        hold on
        ix = year>2001;
        g = (year(ix)-2002)*365+doy(ix);
        plot(xv,xf*splitapply(@mean,fpsn(i,ix),g),'.','Color',cc(i)*ones(1,3))
        xlim([0,730])
        ylim([0,10])
        if i<3
            set(gca,'xtick',cumsum([0,eomday(2001,1:12),eomday(2001,1:12)]))
            set(gca,'xticklabel',years)
        else
            set(gca,'xtick',cumsum([0,eomday(2001,1:12),eomday(2001,1:12)]))
            set(gca,'xticklabel',mls)
            xlabel('Time (month)')
        end
        if i==2||i==4
            set(gca,'yticklabel',[])
        end
        text(10,2,ss((1:3)+(i-1)*3),'FontSize',14,'FontWeight','bold')
        grid on
        box off
    end
    
    ax1 = axes('Position',[0 0 1 0.5],'Visible','off');
    text(0.04,0.39,'GPP (gC/m2/d)',...
        'FontSize',11,'Rotation',90);
    ax1 = axes('Position',[0 0.5 1 0.5],'Visible','off');
    text(0.04,0.2,'Stress Factor (midday)',...
        'FontSize',11,'Rotation',90);
    
    
    xdk.Units = 'Inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    print(xdk,'../figs3/gpp','-dpdf')
    
end

if hh(8)>0
   ix = fsds>400&fsds<425;
   ot = 1:length(fsds);
   xdk = figure;
   cc = get(gca,'ColorOrder');
   co = [2,3,1];
   
          xx = [0.08,0.53];
   
   sssms  = zr*smp(41:60,:);
   sssms  = repmat(sssms(1,mcsec==diurn(10)),48,1);
   sssms  = sssms(1:52551);
   sssms2 = zr*smp(61:80,:);
   sssms2 = repmat(sssms2(1,mcsec==diurn(10)),48,1);
   sssms2 = sssms2(1:52551);
   sssms  = [sssms;sssms2];
ss = {'(a)','(b)'};
   for i = 1:2

       
       ix = fsds>400&fsds<425&year>=2000+i;
             subplot('Position',[xx(i),0.53,0.4,0.41])
       qq = quantile(sssms(i,ix),[1/3,2/3]);
       qq = [-inf,qq,0];
       for j=[3,1,2]
           ix2 = ix&sssms(i,:)>=qq(j)&sssms(i,:)<qq(j+1);
           plot(vpd(ix2),btran(i+2,ix2),'.')
           ylim([0,1])
           hold on
       end
box off
       
       if i==1
           title('AMB')
           ylabel('Stress Factor (SMS)')
       else
           title('TFE')
           set(gca,'yticklabel',[])
       end
       set(gca,'xticklabel',[])
       text(3.5,0.1,ss{i},'FontSize',14,'FontWeight','bold')
   end
   

   
   ssphs  = repmat(vegwp(4,mcsec==diurn(10)),48,1);
   ssphs  = ssphs(1:52551);
   ssphs2 = repmat(vegwp(8,mcsec==diurn(10)),48,1);
   ssphs2 = ssphs2(1:52551);
   ssphs  = [ssphs;ssphs2];
   
   ss = {'(c)','(d)'};
   for i = 1:2
       subplot('Position',[xx(i),0.08,0.4,0.41])
       ix = fsds>400&fsds<425&year>=2000+i;
       qq = quantile(ssphs(i,ix),[1/3,2/3]);
       qq = [-inf,qq,0];
       for j=[3,1,2]
           ix2 = ix&ssphs(i,:)>=qq(j)&ssphs(i,:)<qq(j+1);
           plot(vpd(ix2),btran(i,ix2),'.')
           ylim([0,1])
           hold on
       end
       box off 
       xlabel('VPD (kPa)')
       if i==2
           legend('wettest','driest','intermediate','location','best')
           set(gca,'yticklabel',[])
       else
           ylabel('Stress Factor (PHS)')
       end
       text(3.5,0.1,ss{i},'FontSize',14,'FontWeight','bold')

   end
    
        
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(8)>1
        print(xdk,'../figs3/vpdstress','-dpdf')
    end
    
end



if hh(7)>0

    
    ix = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28); %midday
    g = (year(ix)-2002)*12+month(ix);
    out = splitapply(@mean,btran(:,ix)',g');
    xdk = figure;
       xx = [0.08,0.53];
    for i=1:2
        subplot('Position',[xx(i),0.53,0.4,0.41])
        hold on
        plot(xv,out(:,i+2),'Color',[0.7,0.7,0.7],'LineWidth',2)
        plot(xv,out(:,i),'k','LineWidth',2)
        set(gca,'xtick',2002:0.5:2004)
        set(gca,'xticklabel',[])
        set(gca,'ytick',0:.25:1)
        set(gca,'yticklabel',{'0','','0.5','','1'})
        ylim([0,1])
        if i==1
            title('AMB')
            ylabel('Stress Factor')
            legend({'SMS','PHS'},'Location','SouthEast')

        else
            title('TFE')
            set(gca,'yticklabel',[])
        end
        box off
    end
    
        ix = year>2001;
    g = (year(ix)-2002)*12+month(ix);
    xf = 86400*1e-6*12; %umol/m2/s to g/m2/d
    out = xf*splitapply(@mean,fpsn(:,ix)',g');
    xv  = 2002+(0.5:24)/12;
    
    for i=1:2
        subplot('Position',[xx(i),0.08,0.4,0.41])
        hold on
        plot(xv,out(:,i+2),'Color',[0.7,0.7,0.7],'LineWidth',2)
        plot(xv,out(:,i),'k','LineWidth',2)
        set(gca,'xtick',2002:0.5:2004)
        set(gca,'xticklabel',{'2002','','2003','','2004'})
        set(gca,'ytick',0:2.5:10)
        set(gca,'yticklabel',{'0','','5','','10'})
        ylim([0,10])
        if i==1
            ylabel('GPP (gC/m^2/d)')
        else
            set(gca,'yticklabel',[])
        end
        xlabel('Year')
        
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(7)>0
        print(xdk,'../figs3/gpp','-dpdf')
    end
    
end




if hh(6)>0
    g = (year-2001)*365+doy;
   out  = 4e-7*1800*splitapply(@sum,fctr',g');
   g = (year-2001)*12+month;
   out2 = 4e-7*86400*splitapply(@mean,fctr',g');
   
   xx = [0.08,0.53];
   yy = [0.07,0.365];
   w = 0.4;
   h = 0.27;
   
   dd = [0,0,0;0.7,0.7,0.7];
   mm = {'(PHS)','(SMS)'};
   xdk = figure;
   for i=1:2
   for j=1:2
       ix = p(:,j)>0;
       c = (i-1)*2+j;
       subplot('Position',[xx(j),yy(i),w,h])
       plot(p(ix,j),out(ix,c),'.','Color',dd(i,:))
       
       rr = corr(p(ix,j),out(ix,c))^2;
       text(4.2,1.1,['R^2= ',num2str(round(rr,3))])
       rmse = sqrt(mean((p(ix,j)-out(ix,c)).^2));
       text(4.2,0.5,['RMSE= ',num2str(round(rmse,3))])
       
       std(out(:,c))
       std(p(ix,j))
       
       xlim([0,6])
       ylim([0,6])
       hold on
       box off
       plot([0,6],[0,6],'k:')
       
       if i==2
           set(gca,'xticklabel',[])
       else
           xlabel('Observed Sap Flux (mm/d)')
       end
       if j==2
           set(gca,'yticklabel',[])
       else
           ylabel(['Model ',mm{i}])
       end
   end
   end
   
   
   for i=1:2
       subplot('Position',[xx(i),0.7,w,h]);
       ix = p(:,i)>0;
       ot = 2001+1/365*(0.5:1095)';
       tv = cumsum([0,eomday(2001,1:12)]);
       tv = 0.5*(tv(1:12)+tv(2:13));
       tv = 365+[tv,365+tv];
       tv = 2001+tv/365;
       
       cc=get(gca,'ColorOrder');
       plot(ot(ix),p(ix,i),'.','Color',cc(2,:))
       hold on
       plot(tv,out2(13:36,i+2),'Color',[0.7,0.7,0.7],'LineWidth',1.5)
       plot(tv,out2(13:36,i),'k-','LineWidth',1.5)
       box off
       ylim([0,6])
       set(gca,'xtick',2002:0.5:2004)
       set(gca,'xticklabel',{'2002','','2003','','2004'})
       
       if i==1
           ylabel('Transpiration (mm/d)')
           title('AMB')
       else
           legend('OBS','SMS','PHS')
           title('TFE')
           set(gca,'yticklabel',[])
       end
       xlabel('Year')
       
   end
   
   
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,7];
    xdk.PaperSize = [7,7];
    xdk.PaperPosition = [0,0,7,7];
    
    
    if hh(6)>1
        print(xdk,'../figs3/T','-dpdf')
    end
end


if hh(5)>0
    
    ix1 = smp(1:4,:)==1;
    ix2 = smp(1:4,:)==1;
    ix3 = smp(1:4,:)==1;
    
    x    = vpd;
    y    = 0*smp(1:4,:);
    ffix = fsds>400&fsds<425;
    
    %phs-on, amb
    ee      = 1;
    y(ee,:) = btran(1,:);
    tmp     = repmat(vegwp(4,mcsec==diurn(10)),48,1);
    z       = tmp(1:length(fsds));
    q       = quantile(z(ffix),[1/3,2/3]);
    
    ix1(ee,:) = z<q(1)&ffix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ffix;
    ix3(ee,:) = z>=q(2)&ffix;
    
    
    %phs-on, tfe
    ee      = 2;
    y(ee,:) = btran(2,:);
    tmp     = repmat(vegwp(8,mcsec==diurn(10)),48,1);
    z       = tmp(1:length(fsds));
    ix      = ffix&year>2001;
    q       = quantile(z(ix),[1/3,2/3]);
    
    q/101972
    
    ix1(ee,:) = z<q(1)&ix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ix;
    ix3(ee,:) = z>=q(2)&ix;
    

    ixo=(vpd>1.4&vpd<1.5&fsds>400&fsds<425&btran(2,:)>0.54&btran(2,:)<0.55);
    zz = z(z>=q(2)&ix);
    zz = sort(zz);
    ot = 1:length(zz);
    zo = z(ixo);
    %ot(zz==zo)
    
    fff = fsds(z>=q(2)&ix);
    fff = sort(fff);
    fo  = fsds(ixo);
    %ot(fff==fo)
    
    %phs-off, amb
    ll = 41:60;ee=3;
    y(ee,:) = btran(3,:);
    z = zr*smp(ll,:);
    q = quantile(z(ffix),[1/3,2/3]);
    
    q/101972
    ix1(ee,:) = z<q(1)&ffix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ffix;
    ix3(ee,:) = z>=q(2)&ffix;
    
    %phs-off, amb
    ll = 61:80;ee=4;
    y(ee,:) = btran(4,:);
    z = zr*smp(ll,:);
    ix = ffix&year>2001;
    q = quantile(z(ix),[1/3,2/3]);
    
    ix1(ee,:) = z<q(1)&ix;
    ix2(ee,:) = z>=q(1)&z<q(2)&ix;
    ix3(ee,:) = z>=q(2)&ix;
    
    q/101972
    
    %plotting
    xdk=figure;
    subplot('position',[0.54, 0.56, 0.43, 0.39])
    ee  = 1;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'xticklabel',[])
    set(gca,'ytick',0:0.25:1)
        set(gca,'yticklabel',[])

    text(0.1,0.15,'(a)','FontSize',14,'FontWeight','bold')
    title('TFE')
    
    

        subplot('position',[0.08, 0.12, 0.43, 0.39])
    ee  = 2;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
        set(gca,'ytick',0:0.25:1)
            set(gca,'yticklabel',[])
    xlabel('VPD (kPa)')
    text(0.1,0.15,'(c)','FontSize',14,'FontWeight','bold')

    
    subplot('position',[0.54, 0.12, 0.43, 0.39])
    ee  = 3;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'ytick',0:0.25:1)
    ylabel('Stress function')
    set(gca,'xticklabel',[])
    title('SMS')
    text(2.9,0.15,'(b)','FontSize',14,'FontWeight','bold')
    
    subplot('position',[0.08, 0.56, 0.43, 0.39])
    ee  = 4;
    hold on
    plot(x(ix3(ee,:)),y(ee,ix3(ee,:)),'.','MarkerSize',9)
    plot(x(ix1(ee,:)),y(ee,ix1(ee,:)),'.','MarkerSize',9)
    plot(x(ix2(ee,:)),y(ee,ix2(ee,:)),'.','MarkerSize',9)
    xlim([0 3.25])
    ylim([0 1])
    set(gca,'ytick',0:0.25:1)
    ylabel('Stress function')
    xlabel('VPD (kPa)')
    text(2.9,0.15,'(d)','FontSize',14,'FontWeight','bold')
    title('AMB')
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    
end


if hh(4)>0
    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    
    xx = [0.1,0.55];
    yy = [0.1,0.55];
    w = 0.42;
    h = 0.42;
    
    rst = {'(c)','(d)','(a)','(b)'};
    
    xdk = figure;
    for i=1:2
    for j=1:2
        c = (i-1)*2+j;
        ix = p(:,j)>0;
        subplot('Position',[xx(j),yy(i),w,h])
        plot(p(ix,j),out(ix,c),'r.')
        hold on
        plot([0,5.5],[0,5.5],'k:')
        xlim([0,6])
        ylim([0,6])
        rr = corr(p(ix,j),out(ix,c))^2;
        rmse = sqrt(mean((p(ix,j)-out(ix,c)).^2));
        text(4,0.85,['R^2=',num2str(round(rr,3))])
        text(4,0.5,['RMSE=',num2str(round(rmse,3))])
        text(0.2,5.5,rst{c},'FontSize',14,'FontWeight','Bold')
        box off
        if j==2
            set(gca,'yticklabel',[])
        end
        
        if i==2
            set(gca,'xticklabel',[])
            if j==1
            title('AMB')
            ylabel('SMS T (mm/d)')
            else
            title('TFE')
            end
        else
            if j==1
                ylabel('PHS T (mm/d)')
            end
            xlabel('OBS T (mm/d)')
        end
        
    end
    end
    
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,7];
    xdk.PaperSize = [7,7];
    xdk.PaperPosition = [0,0,7,7];
    
    if hh(4)>0
        print(xdk,'../figs3/trmse','-dpdf')
    end
    
end


if hh(3)>0
   g = (year-2001)*365+doy;
   out = 1800*4e-7*splitapply(@sum,fctr',g');
   
   oneto = (1:1095)';
   ix = oneto>365;
   xv = (2002+(0.5:730)/365)';
   rr = {'(c)','(d)','(a)','(b)','(e)','(f)'};
   
   yy = [0.38,0.69,0.07];
   xx = [0.06,0.53];
   
   xdk = figure;
   ss = [1,2,3,4];
   ii = [1,1,2,2,3,3];
   jj = [1,2,1,2,1,2];
   for i=1:4
       subplot('Position',[xx(jj(i)),yy(ii(i)),0.445,0.28])
       plot(xv,out(ix,i),'k.')
       ylim([0,6])
       set(gca,'xtick',2002:0.5:2004)
       set(gca,'xticklabel',[])
       if i==1

           ylabel('PHS T (mm/d)')
       elseif i==4
           title('TFE')
       elseif i==3
                      title('AMB')
           ylabel('SMS T (mm/d)')
       end
       if i==2||i==4
           set(gca,'yticklabel',[])
       end
       text(2002.05,5.4,rr{i},'FontSize',14,'FontWeight','bold')
   end
   
   
   for i=5:6
       subplot('Position',[xx(jj(i)),yy(ii(i)),0.445,0.28])
       ix = p(:,i-4)>0;
       p(ix,i-4);
       plot(2001+(oneto(ix)-0.5)/365,p(ix,i-4),'k.')
              set(gca,'xtick',2002:0.5:2004)
              if i==5
                  set(gca,'xticklabel',{'2002','','2003','','2004'})
              else
       set(gca,'xticklabel',{'','','2003','','2004'})
              end
       xlabel('Year')
       ylim([0,6])
       text(2002.05,5.4,rr{i},'FontSize',14,'FontWeight','bold')
       if i==5
           ylabel('OBS T (mm/d)')
       else
           set(gca,'yticklabel',[])
       end
   end
   
   
   xdk.Units = 'inches';
   xdk.Position = [2,2,7,6];
   xdk.PaperSize = [7,6];
   xdk.PaperPosition = [0,0,7,6];
   
   if hh(3)>1
        print(xdk,'../figs3/tts','-dpdf')
   end
   
   
end


if hh(2)>0
     
    out = zeros(48,4);
    
    g = findgroups(mcsec);
    for ee=1:4
        ix = (month==9|month==10|month==11) & year==2003;
        out(:,ee) = splitapply(@mean,btran(ee,ix),g(ix));
    end
    
    xdk = figure;
    
    subplot('Position',[0.54,0.11,0.45,0.84])
    hold on
    plot(0.25:0.5:24,out(:,1),'k-','LineWidth',2)
    plot(0.25:0.5:24,out(:,2),'k:','LineWidth',2)
    set(gca,'ytick',0:0.25:1)
    set(gca,'xtick',0:6:24)
    set(gca,'yticklabel',[])
    xlabel('Hour')
    
    legend('AMB','TFE','Location','SouthEast')
    text(1.2,0.1,'(b)','FontSize',14,'FontWeight','bold')
    ylim([0 1])
    title('PHS')
    
    subplot('Position',[0.08,0.11,0.45,0.84])
    
    hold on
    plot(0.25:0.5:24,out(:,3),'k-','LineWidth',2)
    plot(0.25:0.5:24,out(:,4),'k:','LineWidth',2)
    set(gca,'xtick',0:6:24)
    set(gca,'ytick',0:0.25:1)
    title('SMS')
    xlabel('Hour')
    ylabel('Stress factor')
    text(22,0.1,'(a)','FontSize',14,'FontWeight','bold')
    ylim([0 1])
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if hh(2)>1
         print(xdk,'../figs3/fig4','-dpdf')
    end
        
    
end


if hh(1)>0
    
    %figure 2, SON2003 diurnal and monthly mean of veg water potential
    
    % calculate monthly mean midday lwp
    tt=25;
    ixd = mcsec>=diurn(tt)&mcsec<=diurn(tt+3);
    
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
    
    subplot('position',[0.1, 0.59, 0.425, 0.36])
    plot(x,xf*out(1:4,:)')
    xlim([0 24])
    ylim([-3 0])
    xlabel('Hour')
    ylabel({'Water Potential';'(MPa)'})
    title('AMB')
    text(1,-2.5,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:6:24)
    
    subplot('position',[0.545, 0.59, 0.425, 0.36])
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
    
    x = 2001+(0.5:36)/12;
    
    subplot('position',[0.1, 0.1, 0.87, 0.36])
    plot(x,lwp(1,:),'k-','LineWidth',2)
    hold on
    plot(x,lwp(2,:),'k:','LineWidth',2)
    ylim([-3,-1])
    %xlim([0,36])
    set(gca,'ytick',-3:0.5:-1)
    set(gca,'yticklabel',{'-3','','-2','','-1'})
    set(gca,'xtick',2001:0.5:2004)
    set(gca,'xticklabel',{'2001','','2002','','2003','','2004'})
    xlabel('Year')
    ylabel({'Midday Sun Leaf Water';'Potential (MPa)'})
    
    %draw an arrow
    plot(2001+[10.5/12,10.5/12],[-2.6,-2.9],'k-','LineWidth',2)
    plot(2001+[10.5/12,10.2/12],[-2.88,-2.8],'k-','LineWidth',2)
    plot(2001+[10.5/12,10.8/12],[-2.88,-2.8],'k-','LineWidth',2)
    
    legend('AMB','TFE','Location','Southwest')
    text(33.8,-1.3,'(c)','FontSize',14,'FontWeight','bold')
    
      %predawn water potentials?
    out         = xf*out;
    
    disp('AMB MIDDAY')
    disp(mean(out(1:3,25:28),2)')
    
    disp('TFE MIDDAY')
    disp(mean(out(5:7,25:28),2)')
    
    disp('AMB PREDAWN ROOT')
    disp(out(4,11))
    
    disp('TFE PREDAWN ROOT')
    disp(out(8,11))
    
    disp('delta')
    disp(out(4,11)-out(8,11))
    
    disp('AMB MD ROOT & drop')
    disp([mean(out(4,25:28)),mean(out(4,25:28))-out(4,11)])
    
    disp('TFE MD ROOT & drop')
    disp([mean(out(8,25:28)),mean(out(8,25:28))-out(8,11)])
    
    disp('delta drop')
    disp(mean(out(4,25:28))-out(4,11)-(mean(out(8,25:28))-out(8,11)))
    

    disp('AMB MD STEM & drop')
    disp([mean(out(3,25:28)),mean(out(4,25:28))-mean(out(3,25:28))])
    
    disp('TFE MD STEM & drop')
    disp([mean(out(7,25:28)),mean(out(8,25:28))-mean(out(7,25:28))])
    
    disp('delta drop')
    disp(mean(out(4,25:28))-mean(out(3,25:28))-(mean(out(8,25:28))-mean(out(7,25:28))))
    
    
    
    disp('AMB MD SUN & drop')
    disp([mean(out(1,25:28)),mean(out(3,25:28))-mean(out(1,25:28))])
    
    disp('TFE MD SUN & drop')
    disp([mean(out(5,25:28)),mean(out(7,25:28))-mean(out(5,25:28))])
    
    disp('delta drop')
    disp(mean(out(3,25:28))-mean(out(1,25:28))-(mean(out(7,25:28))-mean(out(5,25:28))))
    
    disp('net leaf drop')
    disp(mean(out(1,25:28))-mean(out(5,25:28)))
    
    
    plc1 = 2^(-(mean(out(3,25:28))/-1.75)^2.95);
    plc2 = 2^(-(mean(out(7,25:28))/-1.75)^2.95);
    
    ix = year==2003&month>8&month<12;
    sum(mean(ksr(1:20,ix),2))
    sum(mean(ksr(21:40,ix),2))
    
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if hh(1)>1
        print(xdk,'../figs3/fig2','-dpdf')
    end
    
    
    
end
















































if gg(10)>0
   ix = year==2001 & month==4;
   out = mean(ksr(1:20,ix),2);
   
   plot(z,log10(out),'.')
    
    
end


if gg(9)>0
    
    g = (year-2001)*12+month;
   out =splitapply(@mean,fpsn',g');
   xv = 0.5:36;
   xdk = figure;
    for i=[1,3]
        subplot(2,2,(i==3)+1)
        plot(xv,out(:,i),'k','LineWidth',2)
        hold on
        plot(xv,out(:,i+1),'Color',[0.7,0.7,0.7],'LineWidth',2)
        ylim([0,10])
                xlim([0,36])
        set(gca,'xtick',0:12:36)
        ylabel('GPP')
        if i==1
            title('PHS')
            legend('amb','tfe','location','southwest')
        else
            title('SMS')
        end
    end
    ix  = mcsec>=diurn(25)&mcsec<=diurn(28);
    out =splitapply(@mean,btran(:,ix)',g(ix)');
        for i=[1,3]
        subplot(2,2,(i==3)+3)
        plot(xv,out(:,i),'k','LineWidth',2)
        hold on
        plot(xv,out(:,i+1),'Color',[0.7,0.7,0.7],'LineWidth',2)
        ylim([0,1])
        xlim([0,36])
        set(gca,'xtick',0:12:36)
        ylabel('Stress Factor')
    end
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,6,4];
    xdk.PaperSize = [6,4];
    xdk.PaperPosition = [0,0,6,4];
    
        fout = ['../goodsim/figs_ens/fig6'];
    mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
    
    if gg(9)>0
        print(xdk,fout,'-dpdf')
        system(mkjpg)
    end
    
end

if gg(8)>0
    
    amb_sm = csvread('../goodsim/control_sm.csv');
    amb_sm(amb_sm==0) = nan;
    tfe_sm = csvread('../goodsim/tfe_sm.csv');
    tfe_sm(tfe_sm==0) = nan;
    
    amb_sm(:,1) = 2001+amb_sm(:,1)/365;
    tfe_sm(:,1) = 2001+tfe_sm(:,1)/365;
    
    
    
    %which soil layer
    depths = [nan,0.15,0.5,1,2,3,4,5];
    zzix   = nan(8,1);
    for i=3:8
        j  = 0;
        go = 1;
        while go
            j = j+1;
            if zs(j)>depths(i)
                zzix(i)=j-1;
                go = 0;
            end
        end
    end
    
    dstr = {'0-0.3m','0.5m','1m','2m','3m','4m','5m'};
    tv = 2001+((year-2001)*365+doy)/365;
    pp  = [1,3,2,4];
    
    yyy = {'Soil Water %','Soil Water %','',''};
    ttt = {'PHSamb','PHStfe','SMSamb','SMStfe'};
    for i=2:8
        xdk = figure;
    for j=1:4
        if i==2
            %0-0.3
            dz = zs(2:end)-zs(1:end-1);
            x  = ([dz(1:4),5/6*dz(5)])/0.3;
            targ = x*h2osoi((1:5)+(j-1)*20,:);
        else
            targ = h2osoi(zzix(i)+(j-1)*20,:);
        end
        subplot(2,2,pp(j))
        plot(tv,100*targ,'Color',[0.7,0.7,0.7])
        hold on
        if j==1||j==3
        plot(amb_sm(:,1),amb_sm(:,i),'rx')
        else
            plot(tfe_sm(:,1),tfe_sm(:,i),'rx')
        end
        ylim([5,45])
        xlim([2001,2004])
        ylabel(yyy{j})
        title(ttt{j})
    end
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5,0.97,['Depth = ',dstr{i-1}],'HorizontalAlignment','Center',...
        'FontSize',14,'FontWeight','bold')
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,6,4];
    xdk.PaperSize = [6,4];
    xdk.PaperPosition = [0,0,6,4];
    
    
    fout = ['../goodsim/figs_ens/sm',num2str(i-1)];
    mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];
    
    if gg(8)>0
        print(xdk,fout,'-dpdf')
        system(mkjpg)
    end
    
    
    end
    
end

if gg(7)>0
    ix  = year>2001;
    
    ww  = 7;
    g   = 1+floor((0:length(fctr(ix))-1)/(48*ww));
    out = 1800*4e-7/ww*splitapply(@sum,fctr(:,ix)',g');
    p = [zeros(365,2);p];
    ix1 = p(:,1)>0;
    ix2 = p(:,2)>0;
    ot = 2001+(1:1095)/365;
    
    xv  = 2002+(-ww/2+(1:length(out))*ww)/365;
    xdk = figure;
    
    x = [1,3,2,4];
    for i=1:4
        subplot(2,2,x(i))

        hold on
        if i==1||i==3
        s = scatter(ot(ix1),p(ix1,1),10,[0.7,0.7,0.7],'filled');
                ylabel('AMB T (mm/d)')
        else
        s = scatter(ot(ix2),p(ix2,2),10,[0.7,0.7,0.7],'filled');
                ylabel('TFE T (mm/d)')
        end
        %alpha(s,0.5)
                plot(xv,out(:,i),'r')
        ylim([0,6])
        if i==4
        legend('Obs','Model')
        end
        
        if i==1
            title('PHS')
        elseif i==3
            title('SMS')
        end
        xlim([2002,2004])
    end
    
    
    xdk.Units = 'inches';
xdk.Position = [2,2,6,4];
xdk.PaperSize = [6,4];
xdk.PaperPosition = [0,0,6,4];


fout = '../goodsim/figs_ens/fig5';
mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];

if gg(7)>1
print(xdk,fout,'-dpdf')
system(mkjpg) 
end
    
    
    
end

if gg(6)>0
    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    p = [zeros(365,2);p];
    ix1 = p(:,1)>0;
    ix2 = p(:,2)>0;
    ot = 2001+(1:1095)/365;
    
    xdk = figure;
    
    x = [1,3,2,4];
    for i=1:4
        
        subplot(2,2,x(i))

        if i==1||i==3
            plot(p(ix1,1),out(ix1,i),'.')
            rr = corr(p(ix1,1),out(ix1,i))^2;
            rmse = sqrt(mean((p(ix1,1)-out(ix1,i)).^2));
            hold on
            plot([0,6],[0,6],'k:')
        else
            plot(p(ix2,2),out(ix2,i),'.')
            rr = corr(p(ix2,2),out(ix2,i))^2;
            rmse = sqrt(mean((p(ix2,2)-out(ix2,i)).^2));
            hold on
            plot([0,6],[0,6],'k:')

        end

                    text(0.2,5.5,['R^2=',num2str(round(rr,3))]);
            text(3.5,0.4,['RMSE=',num2str(round(rmse,3))]);
        ylim([0,6])

        
        if i==1
            title('PHS')
            ylabel('Model')
        elseif i==2
            ylabel('Model')
            xlabel('OBS')
        elseif i==3
            title('SMS')
            ylabel('AMB')
        else
            ylabel('TFE')
            xlabel('OBS')
        end
    end
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,6,4];
    xdk.PaperSize = [6,4];
    xdk.PaperPosition = [0,0,6,4];


fout = '../goodsim/figs_ens/fig4';
mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];

if gg(6)>0
print(xdk,fout,'-dpdf')
system(mkjpg)
end    

end

if gg(1)>0
    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    p = [zeros(365,2);p];
    ix1 = p(:,1)>0;
    ix2 = p(:,2)>0;
    ot = 2001+(1:1095)/365;
    
    xdk = figure;
    
    x = [1,3,2,4];
    for i=1:4
        subplot(2,2,x(i))
        plot(ot(ot>=2002),out(ot>=2002,i))
        hold on
        if i==1||i==3
        s = scatter(ot(ix1),p(ix1,1),10,[0.7,0.7,0.7],'filled');
                ylabel('AMB T (mm/d)')
        else
        s = scatter(ot(ix2),p(ix2,2),10,[0.7,0.7,0.7],'filled');
                ylabel('TFE T (mm/d)')
        end
        alpha(s,0.5)
        ylim([0,6])
        if i==2
        legend('Model','Obs')
        end
        
        if i==1
            title('PHS')
        elseif i==3
            title('SMS')
        end
    end
    
    
    xdk.Units = 'inches';
xdk.Position = [2,2,6,4];
xdk.PaperSize = [6,4];
xdk.PaperPosition = [0,0,6,4];


fout = '../goodsim/figs_ens/fig3';
mkjpg = ['convert -density 300 ',fout,'.pdf -quality 100 ',fout,'.jpg'];

print(xdk,fout,'-dpdf')
system(mkjpg)
    
end

if gg(2)>0
    g = (year-2001)*365+doy;
    out = splitapply(@mean,fpsn',g');
    
    for i=1:4
        subplot(2,2,i)
    plot(out(:,i),'.')
    ylim([0,10])
    xlim([0,1095])
    end
    
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    out = splitapply(@mean,btran(:,ix)',g(ix)');
    figure
    for i=1:4
        subplot(2,2,i)
        plot(out(:,i),'.')
        ylim([0,1])
        xlim([0,1095])
    end
    
end

if gg(3)>0
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
    subplot('Position',[0.56,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(a)'};
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
    text(1,-0.5,'(d)','FontWeight','bold','FontSize',14)
    box off
    
    %other 3
    ix = year==2003&month>8&month<12;
    n  = sum(ix);
    b = [0.39,0.1,0.68];
    ct = 0;
    
    for i=1:3
        subplot('Position',[0.08,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end
            if i==2
                xlabel('Day of 2003')
            else
                set(gca,'xticklabel',[])
            end
        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
            set(gca,'xticklabel',[])
        end
        %text(242,0.09*th(i),t{i},'FontWeight','bold')
        text(242,0.09*th(i),t{i})
        text(330,0.01*th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.03,0.39,'Cumulative Water (cm)',...
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

if gg(4)>0
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
    subplot('Position',[0.56,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(a)'};
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
    text(20.5,-0.5,'(d)','FontWeight','bold','FontSize',14)
    box off
    ix = year==2003&month>1&month<5;
    n  = sum(ix);
    b = [0.39,0.1,0.68];
    ct = 0;
    for i=1:3
        subplot('Position',[0.08,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end

        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
        end
        xlim([31,121])
        text(34,th(i),t{i})
        text(90,th(i),p{i},'FontWeight','bold','FontSize',14)
                if i==2
            xlabel('Day of 2003')
        else
                        set(gca,'xticklabel',[])
        end
    end
    ylim([0,110])
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.03,0.39,'Cumulative Water (cm)',...
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


if gg(101)>0

    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    p = [zeros(365,2);p];
    
    ix1 = p(:,1)>0;
    ix2 = p(:,2)>0;
    
    for i=1:4
        subplot(2,2,i)
        if i==1||i==3
            plot(p(ix1,1),out(ix1,i),'.')
            rr = corr(p(ix1,1),out(ix1,i))^2;
            rmse = sqrt(mean((p(ix1,1)-out(ix1,i)).^2));
        else
            plot(p(ix2,2),out(ix2,i),'.')
            rr = corr(p(ix2,2),out(ix2,i))^2;
            rmse = sqrt(mean((p(ix2,2)-out(ix2,i)).^2));
        end
        
        text(4,1.1,['R^2=',num2str(round(rr,3))])
        text(4,0.6,['RMSE=',num2str(round(rmse,3))])
        
        xlim([0,6])
        ylim([0,6])
        hold on
        plot([0,6],[0,6],'k:')
    end
end






ff = [0,0,0,0,0,...
    0,0,0,0,0,...
    0,0,0,0,0,...
    0];


if ff(1)>0
    g = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    
    amb = [zeros(365,1)-1;p(:,1)];
    ix1 = amb>0;
    tfe = [zeros(365,1)-1;p(:,2)];
    ix2 = tfe>0;
    
    subplot(1,2,1)
    plot(out(ix1,1),amb(ix1),'.')
        subplot(1,2,2)
    plot(out(ix2,1),tfe(ix2),'.')
    
    
    corr(out(ix1,1),amb(ix1))^2
    corr(out(ix2,2),tfe(ix2))^2
    sqrt(mean((out(ix1,1)-amb(ix1)).^2))
    sqrt(mean((out(ix2,2)-tfe(ix2)).^2))
end

if ff(2)>0
    ix = year==2003&month>8&month<12;
    g = findgroups(mcsec(ix));
    subplot(1,3,1)
    plot(0.25:0.5:24,splitapply(@mean,smp(23,ix)-vegwp(8,ix)-90,g))
    xlim([6,18])
    grid on
    
    subplot(1,3,2)
    plot(0.25:0.5:24,splitapply(@mean,kon(23,ix),g))
    ylim([0,6e-9])
    xlim([6,18])
    grid on
    
    subplot(1,3,3)
    plot(0.25:0.5:24,splitapply(@mean,qrootsink(23,ix),g))
    xlim([6,18])
    grid on
    
end
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
    subplot('Position',[0.56,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(a)'};
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
    text(1,-0.5,'(d)','FontWeight','bold','FontSize',14)
    box off
    
    %other 3
    ix = year==2003&month>8&month<12;
    n  = sum(ix);
    b = [0.39,0.1,0.68];
    ct = 0;
    
    for i=1:3
        subplot('Position',[0.08,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end
            if i==2
                            xlabel('Day of 2003')
            else
                            set(gca,'xticklabel',[])
            end
        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
                        set(gca,'xticklabel',[])
        end
        %text(242,0.09*th(i),t{i},'FontWeight','bold')
        text(242,0.09*th(i),t{i})
        text(330,0.01*th(i),p{i},'FontWeight','bold','FontSize',14)
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.03,0.39,'Cumulative Water (cm)',...
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
    subplot('Position',[0.56,0.1,0.42,0.83])
    s = {'-',':','-',':'};
    c = [zeros(2,3);0.5*ones(2,3)];
    t = {'Extracted from 0-0.2m','Extracted from 0.2-8.6m','Ambient Precipitation'};
    p = {'(b)','(c)','(a)'};
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
    text(20.5,-0.5,'(d)','FontWeight','bold','FontSize',14)
    box off
    ix = year==2003&month>1&month<5;
    n  = sum(ix);
    b = [0.39,0.1,0.68];
    ct = 0;
    for i=1:3
        subplot('Position',[0.08,b(i),0.42,0.25])
        hold on
        if i<3
            for j=1:4
                ct = ct+1;
                plot(doy(ix)+mcsec(ix)/max(diurn),out2(ct,1:n),...
                    'LineStyle',s{j},'Color',c(j,:),'LineWidth',2)
            end

        else
            plot(doy(ix)+mcsec(ix)/max(diurn),180*cumsum(prec(ix)),'LineWidth',2)
        end
        xlim([31,121])
        text(34,th(i),t{i})
        text(90,th(i),p{i},'FontWeight','bold','FontSize',14)
                if i==2
            xlabel('Day of 2003')
        else
                        set(gca,'xticklabel',[])
        end
    end
    ylim([0,110])
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    text(0.03,0.39,'Cumulative Water (cm)',...
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

if ff(10)>0
    
    %resample to height coords
    out = zeros(860,length(fsds));
    dzn = round(100*(zs(2:end)-zs(1:end-1)));
    ct  = 0;
    for ss=1:20
        out(ct+(1:dzn(ss)),:) = repmat(smp(60+ss,:),dzn(ss),1);
        ct = ct+dzn(ss);
    end
    
    %daily mean
    out = splitapply(@mean,out',(doy+(year-2001)*365)')';
    contour(1:1095,-1:-1:-860,-out/101972,0.1:0.1:2.5)
    
end

if ff(11)>0
    
    ymin = [-1,-1,-3,-3];
    
    figure
    for i=1:4
        t = year+(doy+mcsec/max(diurn)-1)/365;
        dz = 1000*(zs(2:end)-zs(1:end-1));
        subplot(2,2,i)
        plot(t,smp(20*(i-1)+7,:)/101972)
        %plot(t,dz(7)*h2osoi(20*(i-1)+7,:))
        hold on
        plot(t,smp(20*(i-1)+8,:)/101972)
        %plot(t,dz(7)*h2osoi(20*(i-1)+7,:))
        ylim([ymin(i),0])
        %ylim([20,80])
    end

    figure
    for i=1:4
    subplot(2,2,i)
    plot(t,1800*cumsum(qrootsink(7+20*(i-1),:)))
    ylim([0,400])
    end
    
    figure
    out = splitapply(@sum,1800*qrootsink',((year-2001)*12+month)')';
    for i=1:4
        subplot(2,2,i)
    bar(out(7+20*(i-1),:))
    end
end

if ff(12)>0
    
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    
    amb_sm = csvread('../goodsim/control_sm.csv');
    amb_sm(amb_sm==0) = nan;
    tfe_sm = csvread('../goodsim/tfe_sm.csv');
    tfe_sm(tfe_sm==0) = nan;
  
    %which soil layer
    depths = [nan,0.15,0.5,1,2,3,4,5];
    zzix   = nan(8,1);
    for i=3:8
        j  = 0;
        go = 1;
        while go
            j = j+1;
            if zs(j)>depths(i)
                zzix(i)=j-1;
                go = 0;
            end
        end
    end
   

    a = amb_sm;
    pts = [];
    g  = (year-2001)*365+doy;
    
    for i=2:8
        ix1 = ~isnan(a(:,i));
        ix2 = a(ix1,1);
        if sum(ix1)>0
            if i==2
                %0-0.3
                dz = zs(2:end)-zs(1:end-1);
                x  = ([dz(1:4),5/6*dz(5)])/0.3;
                targ = x*h2osoi(1:5,:);
            else
                targ = h2osoi(zzix(i),:);
            end
            out = 100*splitapply(@mean,targ,g);
            pts = [pts;[a(ix1,i),out(ix2)']];
            plot(a(ix1,i),out(ix2),'x')
            
            hold on
        end
    end
    rmse = sqrt(mean((pts(:,1)-pts(:,2)).^2))
    
    xlim([5,45])
    ylim([5,45])
    plot([5,45],[5,45],'k:')
end


if ff(13)>0
    
    tv = doy+mcsec/max(diurn)+(year-2001)*365-1;
    
    a = csvread('../goodsim/control_sm.csv');
    a(a==0) = nan;
    b = csvread('../goodsim/tfe_sm.csv');
    b(b==0) = nan;
    
    p = {'PHSamb','PHStfe','SMSamb','SMStfe'};
    
    dz = -zs(1:5)+zs(2:6);
    dz(end) = dz(end)+0.3-sum(dz);
    
    out = zeros(4,length(fsds));
    for i=1:4
        out(i,:) = dz*h2osoi((1:5)+(i-1)*20,:)/0.3;
    end
    
    
    
    ss = 2;
    dd = 2;
    xdk = figure;
    for i=1:4
        subplot(2,2,i)
        
        if i==1||i==3
            targ = a;
        else
            targ = b;
        end

        
        plot(tv,100*out(i,:))
        hold on
        plot(targ(:,1),targ(:,dd),'rx')
        ylim([5,44])
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
    text(0.5,0.97,'Depth = 0 to 0.3m',...
        'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if ff(13)>0
        print(xdk,'../figs3/soilwater0m','-dpdf')
    end
    
    
    
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

