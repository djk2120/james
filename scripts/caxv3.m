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
        'FCEV','FSH','KSR','ELAI','FGEV','H2OSOI','QOVER','QDRAI'};
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

hh(1:5)   = [0,0,0,0,0];
hh(6:10)  = [0,0,0,0,0];
hh(11:15) = [0,0,0,0,0];
hh(16:20) = [0,0,0,0,0];
hh(21:25) = [0,0,0,0,0];
hh(26:30) = [0,0,0,0,0];
hh(31:35) = [0,0,0,0,0];
hh(36:40) = [0,0,0,0,0];
hh(41:45) = [0,0,0,0,0];
hh(46:50) = [0,0,0,0,0];
hh(51:55) = [0,0,0,0,0];
hh(56:60) = [0,0,0,1,0];



%1  = vwp
%6  = transpiration
%8  = vpd stress
%9  = gpp
%14 = k
%16 = qdry
%17 = qwet
%18 = hr
%19 = smp full profile time series
%20 = h2osoi with obs


if hh(59)>0
    ix = year==2003&month==11;
    for ee = 1:4
        subplot(2,2,ee)
    plot(mean(smp((1:14)+(ee-1)*20,ix),2)/101972,-z(1:14))
    if ee<3
        xlim([-1,0])
    else
        xlim([-3,0])
    end
    end
    
    figure
    ix = year==2003&month>4;
    for ee=1:4
        subplot(2,2,ee)
        plot(splitapply(@mean,smp(7+(ee-1)*20,ix),findgroups(doy(ix)))/101972)
        if ee<3
            ylim([-1,0])
        else
            ylim([-3,0])
        end
    end
    
    figure
    ix = year==2003&month>5;
        for ee=1:4
        subplot(2,2,ee)
        plot(-1+doy(ix)+mcsec(ix)/max(diurn),180*cumsum(qrootsink(7+(ee-1)*20,ix)))
        ylim([0,6])
        xlim([150,365])
        end
    
    
    
    
end

if hh(58)>0
    
    ix = year>2000&month>8&month<12;
    mean(zr*smp(41:60,ix))/101972
    mean(vegwp(4,ix&mcsec==diurn(10)))/101972
    ix = year>2001&month>8&month<12;
    mean(zr*smp(61:80,ix))/101972
    mean(vegwp(8,ix&mcsec==diurn(10)))/101972
    
    
    
    
end


if hh(57)>0
    amb_sm = csvread('../goodsim/control_sm.csv');
    amb_sm(amb_sm==0) = nan;
    tfe_sm = csvread('../goodsim/tfe_sm.csv');
    tfe_sm(tfe_sm==0) = nan;
    
    g  = (year-2001)*365+doy;
    out = splitapply(@mean,h2osoi',g');

    
    yv = -[0.15,0.5,1,2];
    for i=1:18
        subplot(3,6,i)
        plot(tfe_sm(i,2:5)/100,yv,'k-')
        hold on
        plot(out(i,21:32),-z(1:12),'r-')
        plot(out(i,61:72),-z(1:12),'b-')
        title(amb_sm(i,1))
        xlim([0,0.45])
        ylim([-2,0])
    end
    

end


if hh(55)>0
    g = (year-2001)*365+doy;
    
    x = splitapply(@mean,fsds,g);
    for i=1:4
        subplot(2,2,i)
        y = splitapply(@mean,fpsn(i,:),g);
        plot(x,y,'.')
        %xlim([0,2])
        ylim([0,10])
    end
   
    
end


if hh(54)>0
    g = (year-2001)*365+doy;
    
   x = splitapply(@mean,zr*smp(41:60,:),g);
   y = splitapply(@mean,fpsn(3,:),g);
   subplot(2,2,1)
   plot(x/101972,y,'.')
   ylim([0,10])
   xlim([-4,0])
   
   x = splitapply(@mean,zr*smp(61:80,:),g);
   y = splitapply(@mean,fpsn(4,:),g);
   subplot(2,2,2)
   plot(x/101972,y,'.')
   ylim([0,10])
   xlim([-4,0])
   
    
   x = vegwp(4,mcsec==diurn(10));
   y = splitapply(@mean,fpsn(1,:),g);
   subplot(2,2,3)
   plot(x/101972,y,'.')
   ylim([0,10])
   xlim([-4,0])
   
   x = vegwp(8,mcsec==diurn(10));
   y = splitapply(@mean,fpsn(2,:),g);
   subplot(2,2,4)
   plot(x/101972,y,'.')
   ylim([0,10])
   xlim([-4,0])
   
    
end



if hh(53)>0

%how much goes to HR?

pos = smp(1:20,:);
neg = pos;
for i=1:40
    pos(i,:)= qrootsink(i,:).*(qrootsink(i,:)>0);
    neg(i,:)= qrootsink(i,:).*(qrootsink(i,:)<0);
end

   
hrp = 1-sum(qrootsink(1:20,:))./sum(pos(1:20,:));

ix = year==2002;
-180*sum(sum(neg(21:40,ix)))
180*sum(sum(qrootsink(21:40,ix)))


end

if hh(52)>0
    xdk = figure;
    ss = 5;
    
    xv = cell(1,12);
    for i=1:12
        if mod(i,2)==1
            xv{1,i} = -3+i*0.25;
            xv{2,i} = -1.2+i*0.1;
        else
            xv{1,i} = '';
            xv{2,i} = '';
        end
    end
    
    
    %defining k,dp (for SMS)
    ot = 1:n;
    dp = smp;
    for i=41:80
        dp(i,:) = [0,min(189000,smp(i,1:end-1)+255000)];   %model uses last timestep smp!
    end
    dp(dp<100) = nan;
    k = qrootsink./dp;
    k(isnan(k))=0;
    
    ot  = 1:n;
    ix  = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28);
    ix2 = ot(ix)-1;
    subplot(2,2,1)
    plot(smp(ss+60,ix2)/101972,log10(k(ss+60,ix)),'.','Color',...
         [0.8500    0.3250    0.0980])
    hold on
    plot(smp(ss+40,ix2)/101972,log10(k(ss+40,ix)),'.','Color',...
        [0    0.4470    0.7410])
    xlim([-2.75,0])
    ylim([-13,-7])
        set(gca,'xtick',-2.75:0.25:0)
    set(gca,'xticklabel',xv(1,:))
    xlabel('Soil Potential (MPa)')
    ylabel({'Log of (implied)';'conductance [log10(s-1)]'})
    title('SMS')
    box off
    text(10.5/11*-2.75,-13+23/25*6,'(a)','FontSize',14,'FontWeight','bold')
    
    subplot(2,2,2)
    plot(smp(ss+20,ix2)/101972,log10(ksr(ss+20,ix)),'.','Color',...
        [0.8500    0.3250    0.0980])
    hold on
    plot(smp(ss,ix2)/101972,log10(ksr(ss,ix)),'.','Color',...
        [0    0.4470    0.7410])
    xlim([-1.1,0])
    ylim([-13,-7])
    set(gca,'xtick',-1.1:0.1:0)
    set(gca,'xticklabel',xv(2,:))
    xlabel('Soil Potential (MPa)')
    ylabel({'Log of (modeled)';'conductance [log10(s-1)]'})
    title('PHS')
    box off
    text(10.5/11*-1.1,-13+23/25*6,'(b)','FontSize',14,'FontWeight','bold')
    
    
    out = zeros(11,1);
    ix  = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28);
    for i=1:4
        if i>2
            bins = [-2.75:0.25:0];
        else
            bins = [-1.1:0.1:0];
        end
        
        for j= 1:length(bins)-1
            ixb    = ix&smp(ss+(i-1)*20,:)/101972>=bins(j)&smp(ss+(i-1)*20,:)/101972<=bins(j+1);
            out(j,i) = sum(ixb);
        end
        
    end
    
        

    
    subplot(2,2,3)
    b=bar(out(:,3:4));
    b(2).FaceColor = [0.8500    0.3250    0.0980];
    b(1).FaceColor = [0    0.4470    0.7410];
    xlabel('Soil Potential (MPa)')
    xlim([0.5,11.5])
    set(gca,'xtick',0.5:11.5)
    set(gca,'xticklabel',xv(1,:))
    ylabel('Number of timesteps')
    text(1,2300,'(c)','FontSize',14,'FontWeight','bold')
    box off
    
    subplot(2,2,4)
    b=bar(out(:,1:2));
    b(2).FaceColor = [0.8500    0.3250    0.0980];
    b(1).FaceColor = [0    0.4470    0.7410];
    xlim([0.5,11.5])
    set(gca,'xtick',0.5:11.5)
    set(gca,'xticklabel',xv(2,:))
    legend('AMB','TFE','Location','Northwest')
    xlabel('Soil Potential (MPa)')
    ylabel('Number of timesteps')
    text(1,2300,'(d)','FontSize',14,'FontWeight','bold')
    box off
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(52)>1
        print(xdk,'../figs3/suppcond','-dpdf')
    end
    

        
    
    
end


if hh(51)>0
    
    xv = cell(1,12);
    for i=1:12
        if mod(i,2)==1
            xv{i} = -3+i*0.25;
        else
            xv{i} = '';
        end
    end
    
    ix = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28);
    ss = 5;
    
    xdk = figure;
    out = zeros(4,1);

    pp = [3,4,1,2];
    tt = {'','','AMB','TFE'};
    
    xx = [0.04,0.53,0.04,0.53];
    yy = [0.1,0.1,0.55,0.55];
    rr = {'(c)','(d)','(a)','(b)'};
    ee = {'PHS','PHS','SMS','SMS'};
    for i=1:4
        if i>2
            bins = [-2.75:0.25:0];
        else
            bins = [-1.1:0.1:0];
        end
        gg = nan(size(fsds));
        for j=1:length(bins)-1
            ixb = ix&smp(ss+(i-1)*20,:)/101972>=bins(j)&smp(ss+(i-1)*20,:)/101972<=bins(j+1);
            out(i,j) = sum(ixb);
            gg(ixb) = j;
        end
        
        subplot('position',[xx(i),yy(i),0.45,0.4])
        b= boxplot(qrootsink(ss+(i-1)*20,:),gg);
        set(b(7,:),'Visible','off');
        set(b(1:4,:),'Visible','off');
        ylim([0,4e-5])
        title(tt{i})
        aa = get(gca,'xlim');
        if i>2
            xlim([aa(2)-11,aa(2)])
            set(gca,'xtick',aa(2)-11:aa(2))
            text(aa(2)-10.75,3.75e-5,rr{i},'FontSize',14,'FontWeight','bold','VerticalAlignment','middle')
            text(aa(2)-9.75,3.75e-5,'SMS','FontSize',11,'VerticalAlignment','middle')
        else
            xlim([aa(2)-11*2.5,aa(2)])
            set(gca,'xtick',aa(2)-11*2.5:2.5:aa(2))
            text(aa(2)-2.5*10.75,3.75e-5,rr{i},'FontSize',14,'FontWeight','bold')
            text(aa(2)-2.5*9.75,3.75e-5,'PHS','FontSize',11,'VerticalAlignment','middle')
        end
        %text(aa(2)-5,3.7e-5,rr{i},'FontSize',14,'FontWeight','bold')
        %text(aa(2)-48,3.7e-5,ee{i},'FontSize',11,'FontWeight','bold')

        set(gca,'xticklabel',xv)
        
        if i>2
            set(gca,'xticklabel',[])
        else
            xlabel('Soil Potential (MPa)')
        end
        if i==2||i==4
            set(gca,'yticklabel',[])
        else
            ylabel('Layer-5 Midday RWU (mm/s)')
        end
        box off
    end
    
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(51)>1
        print(xdk,'../figs3/rwu','-dpdf')
    end

end

if hh(50)>0
    
    tfe_sm = csvread('../goodsim/tfe_sm.csv');
    tfe_sm(tfe_sm==0) = nan;
    
    subplot(2,2,1)
    plot(100*dz(14)*h2osoi(74,:))
    %ylim([0.05,0.35])
    
        subplot(2,2,2)
    plot(100*dz(14)*h2osoi(34,:))
    %ylim([0.05,0.35])
    
    ix = month>0;
    
    subplot(2,2,3)
    plot(180*cumsum(qrootsink(74,ix)))
    ylim([-1,6])
    
    subplot(2,2,4)
    plot(180*cumsum(qrootsink(34,ix)))
    ylim([-1,6])
    
    figure
    plot(2001+tfe_sm(:,1)/365,tfe_sm(:,6)*dz(14),'rx')
    grid on
    
    
    
end


if hh(49)>0
    ix = year==2003&month>8&month<12;

    yy = [0,6];
    
    subplot(3,2,1)
    bar(180*sum(qrootsink(1:20,ix),2))
    ylim(yy)
    xlim([0,17.5])
    
    subplot(3,2,3)
    bar(180*sum(qrootsink(41:60,ix),2))
    ylim(yy)
    xlim([0,17.5])
    
    subplot(3,2,2)
    bar(180*sum(qrootsink(21:40,ix),2))
    ylim(yy)
    xlim([0,17.5])
    
    subplot(3,2,4)
    bar(180*sum(qrootsink(61:80,ix),2))
    ylim(yy)
    xlim([0,17.5])
    
    subplot(3,2,5)
    bar(zr)
    %ylim([0,6])
    xlim([0,17.5])
    
    subplot(3,2,6)
    bar(zr)
    %ylim([0,6])
    xlim([0,17.5])
    
    180*sum(qrootsink(22,ix))
    figure
    plot(splitapply(@mean,qrootsink(22,ix),findgroups(mcsec(ix))))
    
    figure
    ix = year==2003&doy>285&doy<295;
    xv = doy(ix)+mcsec(ix)/max(diurn);
    plot(xv,qrootsink(2,ix),'.')
    
    %plot(xv,qrootsink(2,ix),'.')
    
end

if hh(48)>0
    xx = [-1,-1,-4,-4];
    g  = (year-2001)*365+doy;
    ss = [vegwp(4,mcsec==diurn(10))/101972;...
        vegwp(8,mcsec==diurn(10))/101972;...
        splitapply(@mean,zr*smp(41:60,:),g)/101972;...
        splitapply(@mean,zr*smp(61:80,:),g)/101972];
    
    ix = p(:,2)>0;
    ll  = fitlm(ss(2,ix),p(ix,2))
    ll2 = fitlm(ss(4,ix),p(ix,2))
    
    
    
end

if hh(47)>0
    oneto = 1:n;
    min(oneto(year==2003&doy==200))
    
    aa = oneto(month==10&year==2003);
    
    for i = min(aa):max(aa)
    plot(smp(41:60,i)/101972,-z,'-x')
    xlim([-3,0])
    ylim([-4,0])
    title(month(i))
    pause(0.01)
    end
    
    
    
    
end


if hh(46)>0
    
    yy=2002;
    xdk = figure;
    n = length(fsds);
    %resample to height
    x = 0;
    ymax = [-1.2,-3];
    vsn  = {'PHS','SMS'};
    pnl  = {'(b)','(a)'};
    cc   = [2,1];
    for ss=[0,40]
        x = x+1;
        nz  = round(100*(zs(2:end)-zs(1:end-1)));
        out = zeros(860,n);
        ct  = 0;
        
        for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) =  repmat(smp(ss+i,:),nz(i),1)/101972;
            ct = ct+nz(i);
        end
        addpath('/Users/kennedy/Documents/MATLAB/othercolor')
        colormap(othercolor('BrBG4'))
        subplot(2,1,cc(x))
        imagesc((1:n)/(48*365)+2001,(1:860)/100,out,[ymax(x),0])
        title(vsn{x})
        xlim([2003,2004])
        ylim([0,4])
        ylabel('Depth (m)')
        if x==1
            xlabel('Month (of 2003)')
        end
        aa={'2003','','','','2004'};
        
        set(gca,'xtick',cumsum([0,eomday(2001,1:12)])/365+2003)
        set(gca,'ytick',0:2:8)
        set(gca,'xticklabel',0:12)
        c = colorbar;
        ylabel(c,'Soil Potential (MPa)','FontSize',11)
        text(2003.02,3.5,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
        
        
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(46)>1
        print(xdk,'../figs3/suppsmp','-djpeg','-r400')
    end
    
end

if hh(45)>0

    yy=2002;
    xdk = figure;
    n = length(fsds);
    %resample to height
    x = 0;
    ymax = [-1.2,-3];
    vsn  = {'PHS','SMS'};
    pnl  = {'(b)','(a)'};
    cc   = [2,1];
    for ss=[20,60]
        x = x+1;
        nz  = round(100*(zs(2:end)-zs(1:end-1)));
        out = zeros(860,n);
        ct  = 0;
        
        for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) =  repmat(smp(ss+i,:),nz(i),1)/101972;
            ct = ct+nz(i);
        end
        addpath('/Users/kennedy/Documents/MATLAB/othercolor')
        colormap(othercolor('BrBG4'))
        subplot(2,1,cc(x))
        imagesc((1:n)/(48*365)+2001,(1:860)/100,out,[ymax(x),0])
        title(vsn{x})
        xlim([2003,2004])
        ylim([0,4])
        ylabel('Depth (m)')
        if x==1
            xlabel('Month (of 2003)')
        end
        aa={'2003','','','','2004'};
        
        set(gca,'xtick',cumsum([0,eomday(2001,1:12)])/365+2003)
        set(gca,'ytick',0:2:8)
        set(gca,'xticklabel',0:12)
        c = colorbar;
        ylabel(c,'Soil Potential (MPa)','FontSize',11)
        text(2003.02,3.5,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
        
        
    end
    
    for i=1:2
        subplot(2,1,i)
        hold on
        plot([2000,2005],[0.5,0.5],'k:')
    end
        
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(45)>1
        print(xdk,'../figs3/smp','-djpeg','-r400')
    end
    
    
    
end


if hh(44)>0
   ix = p(:,1)>0;
   ix2 = p(:,2)>0;
   
   
   
   ixg = mcsec>=diurn(25)&mcsec<=diurn(28);
   g  = (year-2001)*365+doy;
   
   tt = 1800*4e-7*splitapply(@sum,fctr',g');
   
   x1 = splitapply(@mean,fsds,g);
   x2 = splitapply(@mean,vpd,g);
   x3 = splitapply(@mean,btran(1,ixg),g(ixg));
   x4 = splitapply(@mean,btran(2,ixg),g(ixg));
   x = [[x1(ix)';x1(ix2)'],...
       [x2(ix)';x2(ix2)'],...
       [x3(ix)';x4(ix2)']];
   y = [tt(ix,1);tt(ix2,2)];
   
   ll = fitlm(x,y);
   

   subplot(1,2,1)
   a = ll.Coefficients.Estimate;
   bt_resid = y-a(1)-a(2)*x(:,1)-a(3)*x(:,2);
   plot(x(:,3),bt_resid,'.')
   hold on
   plot([0,1],a(4)*[0,1])
   
   x3 = splitapply(@mean,btran(3,ixg),g(ixg));
   x4 = splitapply(@mean,btran(4,ixg),g(ixg));
   x = [[x1(ix)';x1(ix2)'],...
       [x2(ix)';x2(ix2)'],...
       [x3(ix)';x4(ix2)']];
   
   ll2 = fitlm(x,y);
   
   subplot(1,2,2)
   a = ll2.Coefficients.Estimate;
   bt_resid = y-a(1)-a(2)*x(:,1)-a(3)*x(:,2);
   plot(x(:,3),bt_resid,'.')
   hold on
   plot([0,1],a(4)*[0,1],'r-')
   
   
   ll3 = fitlm(x(:,[1,3]),y)
    
    
end



if hh(43)>0
   ix = p(:,1)>0;
   ix2 = p(:,2)>0;
   
   
   
   ixg = mcsec>=diurn(25)&mcsec<=diurn(28);
   g  = (year-2001)*365+doy;
   
  
   
   x1 = splitapply(@mean,fsds,g);
   x2 = splitapply(@mean,vpd,g);
   x3 = splitapply(@mean,btran(1,ixg),g(ixg));
   x4 = splitapply(@mean,btran(2,ixg),g(ixg));
   x = [[x1(ix)';x1(ix2)'],...
       [x2(ix)';x2(ix2)'],...
       [x3(ix)';x4(ix2)']];
   y = [p(ix,1);p(ix2,2)];
   
   ll = fitlm(x,y);
   

   subplot(1,2,1)
   a = ll.Coefficients.Estimate;
   bt_resid = y-a(1)-a(2)*x(:,1)-a(3)*x(:,2);
   plot(x(:,3),bt_resid,'.')
   hold on
   plot([0,1],a(4)*[0,1],'r-')
   
   
   x3 = splitapply(@mean,btran(3,ixg),g(ixg));
   x4 = splitapply(@mean,btran(4,ixg),g(ixg));
   x = [[x1(ix)';x1(ix2)'],...
       [x2(ix)';x2(ix2)'],...
       [x3(ix)';x4(ix2)']];
   
   ll2 = fitlm(x,y);
   
   subplot(1,2,2)
   a = ll2.Coefficients.Estimate;
   bt_resid = y-a(1)-a(2)*x(:,1)-a(3)*x(:,2);
   plot(x(:,3),bt_resid,'.')
   hold on
   plot([0,1],a(4)*[0,1],'r-')
   
    
    
end


if hh(42)>0

    td = csvread('../goodsim/ens_output.csv');
    
    g = (year-2001)*365+doy;
    tt = 1800*4e-7*splitapply(@sum,fctr',g');
    
    sms = zeros(4,1);
    ix = p(:,1)>0;
    sms(1) = sqrt(mean((tt(ix,3)-p(ix,1)).^2));
    sms(2) = corr(tt(ix,3),p(ix,1))^2;
    ix = p(:,2)>0;
    sms(3) = sqrt(mean((tt(ix,4)-p(ix,2)).^2));
    sms(4) = corr(tt(ix,4),p(ix,2))^2;
    
    xdk = figure;
    subplot(1,2,1)
    plot(td(:,8),td(:,9),'.','Color',0.7*ones(1,3))
    hold on
    plot(td(374,8),td(374,9),'rx','MarkerSize',9)
    plot(sms(2),sms(1),'bx','MarkerSize',9)
    ylim([0,2])
    xlim([0,1])
    box off
    ylabel('RMSE (mm/d)')
    xlabel('R^2')
    title('AMB')
    text(0.04,1.92,'(a)','FontSize',14,'FontWeight','bold')
    subplot(1,2,2)
    plot(td(:,12),td(:,13),'.','Color',0.7*ones(1,3))
    hold on
    plot(td(374,12),td(374,13),'rx','MarkerSize',9)
    plot(sms(4),sms(3),'bx','MarkerSize',9)
    ylim([0,2])
    xlim([0,1])
    ylabel('RMSE (mm/d)')
    xlabel('R^2')
    title('TFE')
    box off

    text(0.04,1.92,'(b)','FontSize',14,'FontWeight','bold')
    legend('PHS ensemble',['PHS simulation',newline,'(from main text)'],...
        ['SMS simulation',newline,'(from main text)'])
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];

    if hh(42)>1
    print(xdk,'../figs3/ens','-dpdf')
    end
end


if hh(41)>0
    

    g   = (year-2001)*365+doy;
    tt = 1800*4e-7*splitapply(@sum,fctr',g');
    tt = [tt,p,p];
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    bb = splitapply(@mean,btran(1:2,ix)',g(ix)');
    bb = [bb,splitapply(@mean,btran(3:4,:)',g')];
    bb = [bb,bb];
    pp = [5,7,1,3,6,8,2,4];
    
    for i=1:8
        subplot(4,2,pp(i))
        
        if i<5
            ix = tt(:,i+4)>0;
        else
            ix = tt(:,i)>0;
        end
        plot(bb(ix,i),tt(ix,i),'.')
        xlim([0,1])
        ylim([0,6])
        set(gca,'ytick',0:3:6)
        
    end
    

  
    
end



if hh(40)>0
        g   = (year-2001)*365+doy;
    tt = 1800*4e-7*splitapply(@sum,fctr',g');
    tt = [tt,p,p];
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    bb = splitapply(@mean,btran(1:2,ix)',g(ix)');
    bb = [bb,splitapply(@mean,btran(3:4,:)',g')];
    bb = [bb,bb];

    
        g   = (year-2001)*365+doy;
    tt = 1800*4e-7*splitapply(@sum,fctr',g');
    tt = [tt,p,p];
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    bb = splitapply(@mean,btran(1:2,ix)',g(ix)');
    bb = [bb,splitapply(@mean,btran(3:4,:)',g')];
    bb = [bb,bb];
    pp = [5,7,1,3,6,8,2,4];
    
    
       xx = [0.08,0.55,0.08,0.55];
   yy = [0.1,0.1,0.54,0.54];
   zz = {'(c)','(d)','(a)','(b)'};

       
    
    
    xdk = figure;
    cc = [3,4,1,2];
    for j = 1:4
        bins = 0:0.05:1;
        
        gg = nan(1095,1);
        for i=1:length(bins)-1
            ixb = bb(:,j)>bins(i)&bb(:,j)<=bins(i+1);
            gg(ixb) = i;
        end
        
        subplot('Position',[xx(j),yy(j),0.44,0.4])
        %subplot(2,2,cc(j))
        
        plot([-30,30],[0,0],'k:')
        hold on
        ix = tt(:,j+4)>0;
        
        aa = boxplot(tt(ix,j)-tt(ix,j+4),gg(ix));
        set(aa(7,:),'Visible','off');
        set(aa(1:4,:),'Visible','off');
        
        ylim([-3,3])
        

        aa = get(gca,'xlim');
        xlim([aa(2)-20,aa(2)])
        set(gca,'xtick',aa(2)-20:5:aa(2))
        set(gca,'xticklabel',0:0.25:1)
        
        if j<3
            xlabel('Stress Factor')
        else
            set(gca,'xticklabel',[])
        end
        
        if j==1
            ylabel({'Transpiration (mm/d)';'PHS-obs'})
        end
        
        if j==3
            title('AMB')
            ylabel({'Transpiration (mm/d)';'SMS-obs'})
        end
        
        if j==4
            title('TFE')
        end

        if j==2||j==4
            set(gca,'yticklabel',[])
        end
        
        text(aa(2)-19.4,2.5,zz{j},'FontSize',14,'FontWeight','bold')
        box off
        
        
        
    end
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    if hh(40)>0
        print(xdk,'../figs3/t_v_bt','-dpdf')
    end
    
    gg = zeros(1000,1);
    
end


if hh(39)>0
    g = (year-2001)*365+doy;
    
    bb = splitapply(@mean,btran',g');
    bb = [bb,bb];
    tt = 1800*4e-7*splitapply(@sum,fctr',g');
    tt = [tt,p,p];
        pp = [5,7,1,3,6,8,2,4];
    
    for i=1:8
        if i<5
            ix = tt(:,i+4)>0;
        else
            ix = tt(:,i)>0;
        end
        
        subplot(4,2,pp(i))
        plot(bb(ix,i),tt(ix,i),'.')
        xlim([0,1])
        ylim([0,6])
        set(gca,'ytick',0:3:6)
        
        
        
    end
    
    
end

if hh(38)>0
    ix = p(:,1)>0;
    g  = (year-2001)*365+doy;
    
    tt = [1800*4e-7*splitapply(@sum,fctr',g'),p,p];
    ss = [vegwp(4,mcsec==diurn(10))/101972;...
        vegwp(8,mcsec==diurn(10))/101972;...
        splitapply(@mean,zr*smp(41:60,:),g)/101972;...
        splitapply(@mean,zr*smp(61:80,:),g)/101972];
    
    ss = [ss;ss];
    
    xdk = figure;
    
        
    xx = [-1,-1,-4,-4,-1,-1,-4,-4];
    pp = [5,7,1,3,6,8,2,4];
    
    yy1 = {'PHS','PHS','SMS','SMS','','','',''};
    yy2 = {'amb','tfe','amb','tfe','','','',''};
    
    rr = '(a)(b)(c)(d)(e)(f)(g)(h)';
    
    for i=1:8
        subplot(4,2,pp(i))
        if i<5
            ix = tt(:,i+4)>0;
        else
        ix = tt(:,i)>0;
        end
        plot(ss(i,ix),tt(ix,i),'.')
        ylim([0,6])
        xlim([xx(i),0])
        set(gca,'ytick',0:3:6)
        
        ylabel([yy1{i},yy2{i}])
        
        if pp(i)==1
            title('Model Transpiration (mm/d)')
        end
        if pp(i)==2
            title('Observed Transpiration (mm/d)')
        end
        if pp(i)>6
            xlabel('Model Soil Potential (MPa)')
        end
        box off
        text(xx(i)*0.97,5,rr((1:3)+(pp(i)-1)*3),'FontSize',14,'FontWeight','bold')
    end
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,9];
    xdk.PaperSize = [7,9];
    xdk.PaperPosition = [0,0,7,9];
    
    
    if hh(38)>0
        print(xdk,'../figs3/suppcool','-dpdf')
    end
    
end


if hh(37)>0
    ss = 20;
    yy = 2003;
    for dd=1:365
        mm = mean(month(doy==dd));
   plot(mean(h2osoi((1:20)+ss,year==yy&doy==dd),2),-z,'-x') 
   hold on
   plot(mean(h2osoi((41:60)+ss,year==yy&doy==dd),2),-z,'-x') 
   xlim([0.05,0.35])
   title([num2str(mm),'-',num2str(yy)])
   hold off
    pause(0.1)
    end
    
end


if hh(36)>0
    yy=2002;
    xdk = figure;
    n = length(fsds);
    %resample to height
    x = 0;
    ymax = [-2.5,-2.5];
    vsn  = {'PHS','SMS'};
    pnl  = {'(a)','(b)'};
    for ss=[0,40]
        x = x+1;
        nz  = round(100*(zs(2:end)-zs(1:end-1)));
        out = zeros(860,n);
        ct  = 0;
        
        for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) =  repmat(h2osoi(ss+i,:),nz(i),1);
            ct = ct+nz(i);
        end
        addpath('/Users/kennedy/Documents/MATLAB/othercolor')
        colormap(othercolor('BrBG4'))
        subplot(2,1,x)
        imagesc((1:n)/(48*365)+2001,(1:860)/100,out,[0.1,0.3])
        title(vsn{x})
        xlim([2001,2004])
        ylim([0,8.605])
        ylabel('Depth (m)')
        if x==2
            xlabel('Year')
        end
        aa={'2001','','2002','','2003','','2004'};
        
        set(gca,'xticklabel',aa)
        set(gca,'ytick',0:2:8)
        c = colorbar;
        ylabel(c,'Soil Potential (MPa)','FontSize',11)
        text(2001.08,1,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
        
end


if hh(35)>0
    
    amb_sm = csvread('../goodsim/control_sm.csv');
    amb_sm(amb_sm==0) = nan;
    tfe_sm = csvread('../goodsim/tfe_sm.csv');
    tfe_sm(tfe_sm==0) = nan;
    
    %amb_sm(:,1) = 2001+amb_sm(:,1)/365;
    %tfe_sm(:,1) = 2001+tfe_sm(:,1)/365;
    amb_sm(:,2:end) = amb_sm(:,2:end)/100;
    tfe_sm(:,2:end) = tfe_sm(:,2:end)/100;
    
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
    
    g = (year-2001)*365+doy;

    
    rms = zeros(7,4);
    for j=1:4
        if j==1||j==3
            targ = amb_sm;
            ix = targ(:,1);
        else
            targ = tfe_sm;
            ix = tfe_sm(:,1);
        end
        c = 0;
            ss = 0;
        for i=2:8
            
            if i==2
                x = splitapply(@mean,1/0.3*[dz(1:4),0.1]*h2osoi((1:5)+(j-1)*20,:),g);
            else
                x = splitapply(@mean,h2osoi(zzix(i)+(j-1)*20,:),g);
            end
            
            c = c+1;
            subplot(4,7,c+(j-1)*7)
            plot([0.1,0.4],[0.1,0.4],'k:')
            hold on
            plot(targ(:,i),x(ix),'x')
            xlim([0.1,0.4])
            ylim([0.1,0.4])
            set(gca,'xticklabel',[])
            set(gca,'yticklabel',[])
            
            rms(c,j) = sqrt(mean((targ(:,i)'-x(ix)).^2));
            
            
            out(ss+(1:length(ix)),j) = x(ix);
            if j==1||j==2
                out(ss+(1:length(ix)),j+4) = targ(:,i);
            end
            ss = ss+length(ix);
            
            
        end
    end
    
    
    
    
    figure
    for i=1:4
        subplot(2,2,i)
        if i==1||i==3
            rmse = sqrt(mean((out(:,5)-out(:,i)).^2))
        else
            ix = ~isnan(out(:,6));
            rmse = sqrt(mean((out(ix,6)-out(ix,i)).^2))
        end
    end
    
end


if hh(34)>0
    ix = month>8&month<12;
    
    subplot(2,2,1)
    plot(mean(smp(1:12,ix),2)/101972,-z(1:12),'-x')
    hold on
    plot(mean(smp(41:52,ix)/101972,2),-z(1:12),'-x')
    subplot(2,2,2)
    plot(mean(h2osoi(1:20,ix),2),-z,'-x')
    hold on
    plot(mean(h2osoi(41:60,ix),2),-z,'-x')
    subplot(2,2,3)
    plot(mean(smp(21:32,ix),2)/101972,-z(1:12),'-x')
    hold on
    plot(mean(smp(61:72,ix)/101972,2),-z(1:12),'-x')
    subplot(2,2,4)
    plot(mean(h2osoi(21:40,ix),2),-z,'-x')
    hold on
    plot(mean(h2osoi(61:80,ix),2),-z,'-x')
    
    
    mean(dz*h2osoi(1:20,ix))/mean(dz*h2osoi(41:60,ix))
    mean(dz*h2osoi(21:40,ix))/mean(dz*h2osoi(61:80,ix))
    
    aa = 180*4e-7*sum(fctr,2);
    bb = 180*4e-7*sum(fcev+fgev+qdrai+qover,2);
    
    ix = month>8&month<12&year>2001;
    disp('TFE')
    disp(['SMS: ',num2str(mean(zr*smp(61:80,ix))/101972),' MPa'])
    disp(['PHS: ',num2str(mean(vegwp(8,ix&mcsec==diurn(10)))/101972),' MPa'])
    
end

if hh(33)>0
    yy=2002;
    xdk = figure;
    n = length(fsds);
    %resample to height
    x = 0;
    ymax = [-2.5,-2.5];
    vsn  = {'PHS','SMS'};
    pnl  = {'(b)','(a)'};
    cc   = [2,1];
    for ss=[0,40]
        x = x+1;
        nz  = round(100*(zs(2:end)-zs(1:end-1)));
        out = zeros(860,n);
        ct  = 0;
        
        for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) =  repmat(smp(ss+i,:),nz(i),1)/101972;
            ct = ct+nz(i);
        end
        addpath('/Users/kennedy/Documents/MATLAB/othercolor')
        colormap(othercolor('BrBG4'))
        subplot(2,1,cc(x))
        imagesc((1:n)/(48*365)+2001,(1:860)/100,out,[ymax(x),0])
        title(vsn{x})
        xlim([2001,2004])
        ylim([0,8.605])
        ylabel('Depth (m)')
        if x==2
            xlabel('Year')
        end
        aa={'2001','','2002','','2003','','2004'};
        
        set(gca,'xticklabel',aa)
        set(gca,'ytick',0:2:8)
        c = colorbar;
        ylabel(c,'Soil Potential (MPa)','FontSize',11)
        text(2001.08,1,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    if hh(33)>1
        print(xdk,'../figs3/suppsmp','-dpdf')
    end
        
end



if hh(32)>0
    g = (year-2001)*365+doy;
    qq = 180*splitapply(@sum,qdrai',g');
    

    
    
    
    
    if 1==2
    for yy=2002:2003
        for dd=1:1:365
            mm= mean(month(doy==dd));
            ix = year==yy&doy==dd;
            plot(mean(h2osoi((1:15)+20,ix),2),-z(1:15))
            hold on
            plot(mean(h2osoi((1:15)+60,ix),2),-z(1:15))
            hold off
            xlim([0.05,0.4])
            title([num2str(mm),'-',num2str(yy)])
            pause(0.1)
        end
    end
    end
    
    if 1==2
        for yy=2001:2003
        for dd=1:5:365
            mm= mean(month(doy==dd));
            ix = year==yy&doy==dd;
            plot(mean(smp((1:15)+20,ix),2)/101972,-z(1:15))
            hold on
            plot(mean(smp((1:15)+60,ix),2)/101972,-z(1:15))
            hold off
            xlim([-2.8,0])
            title([num2str(mm),'-',num2str(yy)])
            pause(0.1)
        end
        end
    end
    
    
    
end


if hh(31)>0
    
    yy=2002;
    ix = year==yy&month==9;
    plot(mean(smp((1:15)+20,ix),2)/101972,-z(1:15))
    hold on
    yy=2003;
    ix = year==yy&month==9;
    plot(mean(smp((1:15)+20,ix),2)/101972,-z(1:15))
    
    
    
    if 1==2
    yy = 2002;
    for dd=1:3:365
        plot(mean(smp(1:15,year==yy&doy==dd),2)/101972,-z(1:15))
        hold on
        plot(mean(smp((1:15)+20,year==yy&doy==dd),2)/101972,-z(1:15))
        hold off
        xlim([-0.75,0])
        title(mean(month(doy==dd)))
        pause(0.2)
    end
    end
end
    
    
    
    



if hh(30)>0
    hr = qrootsink(1:40,:);
    for ss=1:40
        hr(ss,:) = -qrootsink(ss,:).*(qrootsink(ss,:)<0);
    end


    180*sum(sum(hr(1:20,year==2003)))
    180*sum(sum(hr(21:40,year==2003)))
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003;
    ix2 = (~(mcsec<diurn(13)|mcsec>diurn(36)))&year==2003;
    
    180*sum(sum(hr(1:20,ix)))
    180*sum(sum(hr(21:40,ix)))
    
end

if hh(29)>0
    
    ix = year==2003&month>8&month<12;
    180*4e-7*sum(fctr(:,ix),2)
    
    if 1==2
    g  = (year-2001)*365+doy;
    xv = 2001+(0.5:1095)/365;
    plot(xv,splitapply(@mean,dz*h2osoi(21:40,:),g)-splitapply(@mean,dz*h2osoi(1:20,:),g))
    hold on
    plot(xv,splitapply(@mean,dz*h2osoi(61:80,:),g)-splitapply(@mean,dz*h2osoi(41:60,:),g))
    end
end


if hh(28)>0
    ix = year==2003&month>8&month<12;
    xv = doy(ix)+mcsec(ix)/max(diurn);
    
    xdk = figure;

    plot(xv,h2osoi(2,ix),'k.')
    hold on
    plot(xv,h2osoi(42,ix),'.','Color',0.5*ones(1,3))
    legend('PHS','SMS')
    ylabel('Volumetric Soil Moisture (-)')
    xlabel('Day of 2003')
    ylim([0.05,0.4])
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if hh(28)>1
        print(xdk,'../figs3/supplayer2','-dpdf')
    end
    
    
end


if hh(27)>0
    dp = smp+255000;
    dp(dp<100)=0;
    k = qrootsink./dp;
    k(dp==0) = nan;
    k(41:60,fctr(:,3)<1) = nan;
    k(61:80,fctr(:,4)<1) = nan;
    k(1:40,:) = ksr(1:40,:);
    
    ot = 1:n;
    ix = mcsec>=diurn(25)&mcsec<=diurn(28)&year==2003;
    xdk = figure;
    rr = {'(a)','(b)','(c)','(d)'};
    tt = {'PHS','SMS'};
    
    i=1;
    subplot(2,2,i)
    plot(smp(25+(i-1)*40,ot(ix)-1)/101972,log10(k(25+(i-1)*40,ot(ix))),'.')
    hold on
    plot(smp(5+(i-1)*40,ot(ix)-1)/101972,log10(k(5+(i-1)*40,ot(ix))),'.')
    ylim([-13,-7])
    xlim([-2.5,0])
    title(tt{i})
    ylabel({'Log-10 of Layer 5';'conductance [log(1/s)]'})
    box off
    text(-2.4,-12.3,rr{i},'FontSize',14,'FontWeight','bold')
    
    i=2;
    subplot(2,2,i)
    plot(smp(25+(i-1)*40,ot(ix))/101972,log10(k(25+(i-1)*40,ot(ix))),'.')
    hold on
    plot(smp(5+(i-1)*40,ot(ix))/101972,log10(k(5+(i-1)*40,ot(ix))),'.')
    ylim([-13,-7])
    xlim([-2.5,0])
    title(tt{i})
    box off
    text(-2.4,-12.3,rr{i},'FontSize',14,'FontWeight','bold')
    
    i=1;
    subplot(2,2,i+2)
    if i==2
        plot([-10,1],[0,0],'k:')
    end
    hold on
    plot(smp(25+(i-1)*40,ot(ix)-1)/101972,qrootsink(25+(i-1)*40,ix),'.')
    
    plot(smp(5+(i-1)*40,ot(ix)-1)/101972,qrootsink(5+(i-1)*40,ix),'.')
    ylim([-2e-5,5e-5])
    box off
    xlim([-2.5,0])
    if i==1
        plot([-10,1],[0,0],'k:')
        ylabel({'Layer 5 Root Water Uptake';'(mm/s)'})
        legend('TFE','AMB','location','northwest')
    end
    text(-2.4,7/60*7e-5-2e-5,rr{i+2},'FontSize',14,'FontWeight','bold')
    
        i=2;
    subplot(2,2,i+2)
    if i==2
        plot([-10,1],[0,0],'k:')
    end
    hold on
    plot(smp(25+(i-1)*40,ix)/101972,qrootsink(25+(i-1)*40,ix),'.')
    
    plot(smp(5+(i-1)*40,ix)/101972,qrootsink(5+(i-1)*40,ix),'.')
    ylim([-2e-5,5e-5])
    box off
    xlim([-2.5,0])
    if i==1
        plot([-10,1],[0,0],'k:')
        ylabel({'Layer 5 Root Water Uptake';'(mm/s)'})
        legend('TFE','AMB','location','northwest')
    end
    text(-2.4,7/60*7e-5-2e-5,rr{i+2},'FontSize',14,'FontWeight','bold')

    
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(27)>1
        print(xdk,'../figs3/suppcond','-dpdf')
    end
    
    
end



if hh(26)>0

   g  = (year-2001)*365+doy; 
    ss = [vegwp(4,mcsec==diurn(10))/101972;...
        vegwp(8,mcsec==diurn(10))/101972;...
        splitapply(@mean,zr*smp(41:60,:),g)/101972;...
        splitapply(@mean,zr*smp(61:80,:),g)/101972];
   
   tt = 1800*4e-7*splitapply(@sum,fctr',g');
   rr = {'AMB','TFE'};
   kk = [1,2,1,2];
   cc = [3,4,1,2];
   
   xx = [0.08,0.55,0.08,0.55];
   yy = [0.1,0.1,0.55,0.55];
   zz = {'(c)','(d)','(a)','(b)'};
   
   xdk = figure;
   for j = 1:4
       k = kk(j);
       subplot('Position',[xx(j),yy(j),0.44,0.38])
       
       ix = p(:,k)>0;
       gg = nan(size(p(:,1)));
       if j<3
           bins = -1:0.05:0;
       else
           bins = -4:0.2:0;
       end
       out = zeros(8,1);
       for i = 1:length(bins)-1
           ixb = ix'&ss(j,:)>bins(i)&ss(j,:)<=bins(i+1);
           gg(ixb) = i;
           out(i) = mean(tt(ixb,j)-p(ixb,k));
       end
       plot([-100,100],[0,0],'k:')
       hold on
       bb= boxplot(tt(:,j)-p(:,k),gg);
       set(bb(7,:),'Visible','off');
       set(bb(1:4,:),'Visible','off');
       ylim([-3,3])
       
       aa = get(gca,'xlim');
       xlim([aa(2)-20,aa(2)])
       
       if j==2||j==4
           set(gca,'yticklabel',[])
       elseif j==1
           ylabel({'Transpiration (mm/d)';'PHS-Obs'})
       else
           ylabel({'Transpiration (mm/d)';'SMS-Obs'})
       end
       
       if j<3
           set(gca,'xtick',aa(2)-20:5:aa(2))
           set(gca,'xticklabel',-1:0.25:0)
           xlabel('Model Soil Potential (MPa)')
       else
           set(gca,'xtick',aa(2)-20:5:aa(2))
           set(gca,'xticklabel',-4:1:0)
           title(rr{j-2})
           
       end
       text(aa(2)-19.4,2.5,zz{j},'FontSize',14,'FontWeight','bold')
       box off
       
   end
   ax1 = axes('Position',[0 0 1 1],'Visible','off');
   %text(0.514,0.076,'0','FontSize',10)
   %text(0.9832,0.076,'0','FontSize',10)

               xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    if hh(26)>0
        print(xdk,'../figs3/sm3','-dpdf')
    end
    
end

if hh(25)>0
    ix  = mcsec>0;
    g   = month(ix);
    out = splitapply(@mean,fpsn(:,ix)',g');
    subplot(1,2,1)
    bar(out(:,1)-out(:,3))
    ylim([-2.5,2.5])
    
    ix  = mcsec>=diurn(25)&mcsec<=diurn(28);
    g   = month(ix);
    [(1:12)',splitapply(@mean,btran(1,ix)',g')]
    
    ix  = year>2001;
    g   = month(ix);
    out = splitapply(@mean,fpsn(:,ix)',g');
    subplot(1,2,2)
    bar(out(:,2)-out(:,4))
    ylim([-2.5,2.5])
    
    ix  = year>2001&mcsec>=diurn(25)&mcsec<=diurn(28);
    g   = month(ix);
    [(1:12)',splitapply(@mean,btran(2,ix)',g')]
    
    ix = year>2001;
    g   = month(ix);
    [(1:12)',splitapply(@mean,vpd(ix)',g')]
    
    ix = year>2001&mcsec==diurn(10);
    g   = month(ix);
    [(1:12)',splitapply(@mean,vegwp(8,ix)',g')/101972]
    
    
end


if hh(24)>0
    g   = (year-2001)*365+doy;
    out = 1800*4e-7*splitapply(@sum,fctr',g');
    
    ss = vegwp(4,mcsec==diurn(10));
    ss = [ss;vegwp(8,mcsec==diurn(10))];
    ss = [ss;splitapply(@mean,zr*smp(41:60,:),g)];
    ss = [ss;splitapply(@mean,zr*smp(61:80,:),g)];
    ss = ss/101972;
    
    xdk = figure;
    pp = [1,2,1,2];
    xx = [-1,-1,-4,-4];
    rr = {'(c)','(d)','(a)','(b)'};
    cc = [3,4,1,2]
    for i=1:4
        subplot(2,2,cc(i))
        ix = p(:,pp(i))>0;
        plot(ss(i,ix),out(ix,i)-p(ix,pp(i)),'k.')
        xlim([xx(i),0])
        ylim([-3,3])
        
        if i==1
            ylabel({'Transpiration (mm/d)';'PHS-Obs'})
        end
        
        if i<3
            xlabel('Model Soil Potential (MPa)')
        end
        
        if i==3
            ylabel({'Transpiration (mm/d)';'SMS-Obs'})
            title('AMB')
        end
        if i==4
            title('TFE')
        end
        if i==2||i==4
            set(gca,'yticklabel',[])
        end
        box off
        text(0.97*xx(i),2.6,rr{i},'FontSize',14,'FontWeight','bold')
    end
    
            xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];

    if hh(24)>0
        print(xdk,'../figs3/t_vs_smp','-dpdf')
    end
    
    
end


if hh(23)>0
    zz = 0:860;
    zr = 0.98.^(zz);
    
    xdk = figure;
    plot(zr,0:-0.01:-8.6,'k','LineWidth',2)
    ylim([-9,0])
    set(gca,'ytick',-9:1:0)
    xx = cell(10,1);
    xx(:) = {''};
    xx(1:3:10) = [{'9'},{'6'},{'3'},{'0'}];
    set(gca,'yticklabel',xx)
    ylabel('Depth (m)')
    xlabel('Cumulative Root Fraction')
    box off
    grid on
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,4,3];
    xdk.PaperSize = [4,3];
    xdk.PaperPosition = [0,0,4,3];

    if hh(23)>0
        print(xdk,'../figs3/roots','-dpdf')
    end
    
end



if hh(22)>0
    n = length(fsds);
    out = zeros(2,n);
    for ee=0:1
        for i=1:length(prec)
            t = qrootsink((1:20)+ee*20,i);
            surp = 0;
            up   = 0;
            down = 0;
            for ss=20:-1:1
                x = t(ss);
                if x>0
                    surp = surp+x;
                elseif ss==2
                    up = up-x;
                elseif abs(x)<=surp
                    up = up-x;
                    surp = surp+x;
                elseif surp>0
                    up = up+surp;
                    down = down-x-surp;
                    surp = surp+x;
                else
                    surp = surp+x;
                    down = down-x;
                    
                end
            end
            out(1+ee*2,i)=up;
            out(2+ee*2,i)=down;
        end
    end
    
    xdk=figure;

    ee = [0,1,0,1];
    yy = [2002,2002,2003,2003];
    tt = {'AMB-2002','TFE-2002','AMB-2003','TFE-2003'};
    
    addpath('/Users/kennedy/Documents/MATLAB/othercolor')
    cc=colormap(othercolor('BrBG4'));
    
    rr = {'(a)','(b)','(c)','(d)'};
    
    for i=1:4
        subplot(2,2,i)
        ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==yy(i);
        ix2 = (~(mcsec<diurn(13)|mcsec>diurn(36)))&year==yy(i);
        a1=(splitapply(@sum,180*out(1+ee(i)*2,ix),month(ix)));
        a2=(splitapply(@sum,180*out(1+ee(i)*2,ix2),month(ix2)));
        b1=(splitapply(@sum,180*out(2+ee(i)*2,ix),month(ix)));
        b2=(splitapply(@sum,180*out(2+ee(i)*2,ix2),month(ix2)));
        
        b=bar(2:3:36,a1+a2,0.3);
        b.FaceColor=cc(20,:);
        b.EdgeColor=cc(15,:);
        hold on
        b=bar(2:3:36,a1,0.3);
        b.FaceColor=cc(5,:);
        b.EdgeColor=cc(4,:);
        b=bar(1:3:36,b1+b2,0.3);
        b.FaceColor=cc(40,:);
        b.EdgeColor=cc(50,:);
        b=bar(1:3:36,b1,0.3);
        b.FaceColor=cc(60,:);
        b.EdgeColor=cc(61,:);
        xlim([0,36])
        ylim([0,10])
        set(gca,'xtick',1.5:3:36)
        set(gca,'xticklabel',1:12)
        xlabel('Month')
        ylabel('HR (cm)')
        box off
        title(tt{i})
        text(0.9,9.3,rr{i},'FontSize',14,'FontWeight','bold')
        if i==2
                legend('Up,day','Up,night',...
        'Down,day','Down,night')
        end
    end

    % out = up,down,up,down
    % out = a,a,tfe,tfe
    
    180*sum(sum(out(1:2,year==2003)))
    180*sum(sum(out(3:4,year==2003)))
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];

    if hh(22)>0
        print(xdk,'../figs3/hr2','-dpdf')
    end
    
end

if hh(21)>0
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
    
    g = (year-2001)*365+doy;
    xv = 2001+(0.5:1095)/365;
    
    
    tt = {'AMB','TFE'};
    xx = [0.12,0.57,0.12,0.57];
    yy = [0.08,0.08,0.53,0.53];
    rr = {'(c)','(d)','(a)','(b)'};

    xdk = figure;
    ss = [2,4,1,3];
    
    uu = [];
    
    for i=[2,4:8]
    for j=1:4

        if i>3
            jj=subplot(6,4,ss(j)+(i-3)*4);
            plot(xv,splitapply(@mean,h2osoi(zzix(i)+(j-1)*20,:),g),'LineWidth',1.5)
            hold on
        else

            jj=subplot(6,4,ss(j));
            plot(xv,splitapply(@mean,1/0.3*[dz(1:4),0.1]*h2osoi((1:5)+(j-1)*20,:),g),'LineWidth',1.5)
            hold on
        end
        if j==1||j==3
            plot(amb_sm(:,1),amb_sm(:,i)/100,'rx')
        else
            plot(tfe_sm(:,1),tfe_sm(:,i)/100,'rx')
        end
        ylim([0.05,0.35])
        xlim([2001,2004])
        
        set(gca,'xtick',2001:2004)

        if i~=8
            set(gca,'xticklabel',[])
        end
        
        
        set(gca,'ytick',0.05:0.1:0.35)
        if j~=3
            set(gca,'yticklabel',[])
        end
        uu = [uu,{jj}];
    
        box off
    end
    end
    
    legend('Model','Obs')
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    
    yy = [0,0.88,0];
    yy = [yy,yy(2)-0.1425:-0.1425:0];
    dd = {'','0-0.3m','','1m','2m','3m','4m','5m'};
    for i=[2,4:8]
    text(0.9,yy(i),dd{i},'FontSize',11,'FontWeight','bold')
    end
    dd = {'SMSamb','PHSamb','SMStfe','PHStfe'};
    xx = [0.205:0.2062:1];
    for i=1:4
        text(xx(i),0.95,dd{i},'FontSize',12,'FontWeight','bold',...
            'HorizontalAlignment','Center')
    end
    text(0.05,0.42,'Volumetric Soil Moisture (-)','FontSize',11,...
        'Rotation',90)
    text(0.5,0.05,'Year','FontSize',11,...
        'HorizontalAlignment','Center')
        
        
    
    
        xdk.Units = 'inches';
    xdk.Position = [2,2,7,9];
    xdk.PaperSize = [7,9];
    xdk.PaperPosition = [0,0,7,9];
    if hh(21)>1
        print(xdk,'../figs3/suppsm2','-dpdf')
    end
    
    
    
end

if hh(20)>0
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
    
    g = (year-2001)*365+doy;
    xv = 2001+(0.5:1095)/365;
    i=3;
    
    tt = {'AMB','TFE'};
    xx = [0.12,0.57,0.12,0.57];
    yy = [0.08,0.08,0.53,0.53];
    rr = {'(c)','(d)','(a)','(b)'};

    xdk = figure;
    
    ss = [3,4,1,2];
    for j=1:4
        subplot('Position',[xx(j),yy(j),0.4,0.41])
        plot(xv,splitapply(@mean,h2osoi(zzix(3)+(j-1)*20,:),g),'LineWidth',1.5)
        hold on
        if j==1||j==3
            plot(amb_sm(:,1),amb_sm(:,3)/100,'rx')
        else
            plot(tfe_sm(:,1),tfe_sm(:,3)/100,'rx')
        end
        ylim([0.05,0.35])
        
        set(gca,'xtick',2001:2004)
        if ss(j)<3
        title(tt{ss(j)})
        set(gca,'xticklabel',[])
        else
            xlabel('Year')
        end
        
        text(2001.1,0.08,rr{j},'FontSize',14,'FontWeight','bold')
        set(gca,'ytick',0.05:0.1:0.35)
        
        if j==1
            ylabel({'PHS (depth = 50cm)';'Volumetric Soil Moisture (-)'})
        elseif j==3
            ylabel({'SMS (depth = 50cm)';'Volumetric Soil Moisture (-)'})
        else
            set(gca,'yticklabel',[])
        end

        box off
        
        if j==2||j==4
                plot(2001+[10/12,10/12],[.3*4/30+.05,.3*1/30+.05],'k-','LineWidth',2)
    plot(2001+[10/12,9.7/12],[.3*12/300+.05,.2*4/30+.05],'k-','LineWidth',2)
    plot(2001+[10/12,10.3/12],[.3*12/300+.05,.2*4/30+.05],'k-','LineWidth',2)
        end
                if j==2
            legend({'Model','Obs'},'Location','Northeast')
        end
    end
    
    
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    if hh(20)>1
        print(xdk,'../figs3/sm2','-dpdf')
    end
    
    
    
end

if hh(19)>0
    yy=2002;
    xdk = figure;
    n = length(fsds);
    %resample to height
    x = 0;
    ymax = [-2.5,-2.5];
    vsn  = {'PHS','SMS'};
    pnl  = {'(b)','(a)'};
    cc   = [2,1];
    for ss=[20,60]
        x = x+1;
        nz  = round(100*(zs(2:end)-zs(1:end-1)));
        out = zeros(860,n);
        ct  = 0;
        
        for i=1:20
            ix = ct+(1:nz(i));
            out(ix,:) =  repmat(smp(ss+i,:),nz(i),1)/101972;
            ct = ct+nz(i);
        end
        addpath('/Users/kennedy/Documents/MATLAB/othercolor')
        colormap(othercolor('BrBG4'))
        subplot(2,1,cc(x))
        imagesc((1:n)/(48*365)+2001,(1:860)/100,out,[ymax(x),0])
        title(vsn{x})
        xlim([2001,2004])
        ylim([0,8.605])
        ylabel('Depth (m)')
        if x==2
            xlabel('Year')
        end
        aa={'2001','','2002','','2003','','2004'};
        
        set(gca,'xticklabel',aa)
        set(gca,'ytick',0:2:8)
        c = colorbar;
        ylabel(c,'Soil Potential (MPa)','FontSize',11)
        text(2001.08,1,pnl{x},'FontSize',14,'FontWeight','bold','Color',[0.8,0.8,0.8])
    end
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    if hh(19)>1
        print(xdk,'../figs3/smp','-dpdf')
    end
        
end


if hh(18)>0
    
    hr = qrootsink(1:40,:);
    for ss=1:40
        hr(ss,:) = -qrootsink(ss,:).*(qrootsink(ss,:)<0);
    end
    
    ix = (mcsec<diurn(13)|mcsec>diurn(36))&year==2003;
    ix2 = (~(mcsec<diurn(13)|mcsec>diurn(36)))&year==2003;
    hr1 = [180*splitapply(@sum,sum(hr(1:20,ix)),month(ix))',...
        180*splitapply(@sum,sum(hr(1:20,ix2)),month(ix2))'];
    hr2 = [180*splitapply(@sum,sum(hr(21:40,ix)),month(ix))',...
        180*splitapply(@sum,sum(hr(21:40,ix2)),month(ix2))'];
    
    
    
    xdk = figure;
    
    subplot('Position',[0.08,0.11,0.45,0.84])
    b=bar(hr1,'stacked');
    b(1).FaceColor = [0.4,0.4,0.6];
    b(2).FaceColor = [0.6,0.6,0.8];
    b(1).EdgeColor = [0.35,0.35,0.6];
    b(2).EdgeColor = [0.55,0.55,0.8];
    xlabel('Month')
    ylabel('Total HR (cm)')
    title('AMB')
    xlim([0 13])
    ylim([0 10])
    text(0.4,9.5,'(a)','FontSize',14,'FontWeight','bold')
    
    
    subplot('Position',[0.54,0.11,0.45,0.84])
    b=bar(hr2,'stacked');
    b(1).FaceColor = [0.4,0.4,0.6];
    b(2).FaceColor = [0.6,0.6,0.8];
    b(1).EdgeColor = [0.35,0.35,0.6];
    b(2).EdgeColor = [0.55,0.55,0.8];
    %bar(hr2,'FaceColor',[0.6,0.6,0.8],'EdgeColor',[0.55,0.55,0.8])
    xlabel('Month')
    xlim([0 13])
    ylim([0 10])
    title('TFE')
    set(gca,'yticklabel',[])
    text(0.4,9.5,'(b)','FontSize',14,'FontWeight','bold')
    legend({'Nighttime','Daytime'})
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if hh(18)>1
        print(xdk,'../figs3/hr','-dpdf')
    end
    
end


if hh(17)>0
    
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
    ymax = [40,20,110];
    ymin = [0,-20,0];
    tlab   = {'PHS-amb','PHS-tfe','SMS-amb','SMS-tfe'};
    for i=1:4
        plot(cumsum(flipud(out(1:400,i))),-4:0.01:-0.01,...
            'LineStyle',s{i},'Color',c(i,:),'LineWidth',2)
        hold on
    end
    xlim([-15,25])
    set(gca,'ytick',-4:1:0)
    set(gca,'yticklabel',4:-1:0)
    ylabel('Depth (m)')
    xlabel({'Cumulative Root Water Uptake (cm)'})
    legend(tlab,'location','Southeast')
    text(2/30*40-15,-3.75,'(d)','FontWeight','bold','FontSize',14)
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
        ylim([ymin(i),ymax(i)])
        text(34,ymin(i)+0.9*(ymax(i)-ymin(i)),t{i})
        text(112,ymin(i)+0.12*(ymax(i)-ymin(i)),p{i},'FontWeight','bold','FontSize',14)
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
    
    if hh(17)>1
        print(xdk,'../figs3/qwet','-dpdf')
    end
    

    go = 1;
    i  = 10;
    while go
        i = i+1;
        if sum(out(i:end,1))<0
            i
            go = 0;
        end
        if i==859
            i
            go = 0;
        end
    end
    sum(out(i-1:end,1))
    sum(out(i:end,1))
    dd = i-1+sum(out(i-1:end,1))/out(i-1,1)
    
    
    disp('SMS extracts:')
    ttt = (sum(out(i:end,3))+out(i-1,3)*sum(out(i-1:end,1))/out(i-1,1));
    disp(ttt)
    disp(ttt/sum(out(:,3)))
    
    
    go = 1;
    i  = 1;
    while go
        i = i+1;
        if sum(out(i:end,2))<0
            i
            go = 0;
        end
        if i==859
            i
            go = 0;
        end
    end
    sum(out(i-1:end,2))
    sum(out(i:end,2))
    sum(out(i-1:end,2))/out(i-1,2)
    dd = i-1+sum(out(i-1:end,2))/out(i-1,2)
    
    disp('SMS extracts:')
    ttt = (sum(out(i:end,4))+out(i-1,4)*sum(out(i-1:end,2))/out(i-1,2));
    disp(ttt)
    disp(ttt/sum(out(:,4)))
        
        
    
end


if hh(16)>0
    
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
    t = {'Extracted from 0-0.2m','Extracted beyond 0.2m','Ambient Precipitation'};
    p = {'(b)','(c)','(a)'};
    tlab   = {'PHS-amb','PHS-tfe','SMS-amb','SMS-tfe'};
    
    ymax = [20,20,30];
    
    for i=1:4
        plot(cumsum(flipud(out(1:400,i))),-4:0.01:-0.01,...
            'LineStyle',s{i},'Color',c(i,:),'LineWidth',2)
        hold on
    end
    set(gca,'ytick',-4:0)
    set(gca,'yticklabel',4:-1:0)
    ylabel('Depth (m)')
    xlabel({'Cumulative Root Water Uptake (cm)'})
    legend(tlab,'location','Southeast')
    text(2,-3.75,'(d)','FontWeight','bold','FontSize',14)
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
        
        ylim([0,ymax(i)])
        text(242,0.9*ymax(i),t{i})
        text(330,0.12*ymax(i),p{i},'FontWeight','bold','FontSize',14)
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
    
    if hh(16)>1
        print(xdk,'../figs3/qdry','-dpdf')
    end
end



if hh(15)>0
    
    dp = smp+255000; %sms
    dp(1:20,:) = smp(1:20,:)-repmat(vegwp(4,:),20,1); %phsamb (ignores gravity)
    dp(21:40,:) = smp(21:40,:)-repmat(vegwp(8,:),20,1); %phstfe (ignores gravity)
    %layer 5 gravity ~ 0.0025 MPa
    
    xdk = figure;
    ix = year==2003;
    g = findgroups(mcsec(ix));
    xv = 0.25:0.5:24;
    yy = [0.3,2.5];
    yy2 = [-0.15,0];
    cc  = [0,0.7];
    tt  = {'PHS','SMS'};
    pp  = {'(b)','(a)'};
    zz  = [2.35*45/250-0.15,2.35];
    for i=1:2
        subplot(1,2,(i==1)+1)
        hold on
        plot(xv,splitapply(@mean,dp(5+(i-1)*40,ix),g)/101972,'-','LineWidth',2,'Color',cc(i)*ones(1,3))
        plot(xv,splitapply(@mean,dp(5+(i-1)*40+20,ix),g)/101972,':','LineWidth',2,'Color',cc(i)*ones(1,3))
        
        xlim([0,24])
        set(gca,'xtick',0:6:24)
        ylim([yy2(i),yy(i)])
        title(tt{i})
        xlabel('Hour of Day')
        ylabel({'Average Layer 5 \Delta\psi (MPa)'})
        text(1.1,zz(i),pp{i},'FontSize',14,'FontWeight','bold')
        if i==2
            legend('AMB','TFE','Location','SouthWest')
        end
        
    end
    
    xdk.Units = 'Inches';
    xdk.Position = [2,2,7,4];
    xdk.PaperSize = [7,4];
    xdk.PaperPosition = [0,0,7,4];
    
    if hh(15)>1
        print(xdk,'../figs3/suppdp','-dpdf')
    end
    
end


if hh(14)>0
    
    %defining k,dp (for SMS)
    ot = 1:n;
    dp = smp;
    for i=41:80
        dp(i,:) = [0,min(189000,smp(i,1:end-1)+255000)];   %model uses last timestep smp!
    end
    dp(dp<100) = nan;
    k = qrootsink./dp;
    k(isnan(k))=0;
    
    %averaging
    ix = year==2003;
    out = [splitapply(@mean,ksr(5,ix),findgroups(mcsec(ix)));...
        splitapply(@mean,ksr(25,ix),findgroups(mcsec(ix)));...
        splitapply(@mean,k(45,ix),findgroups(mcsec(ix)));...
        splitapply(@mean,k(65,ix),findgroups(mcsec(ix)))];
    out2 = [splitapply(@mean,smp(5,ix)-vegwp(4,ix),findgroups(mcsec(ix)));...
        splitapply(@mean,smp(25,ix)-vegwp(8,ix),findgroups(mcsec(ix)));...
        splitapply(@mean,min(189000,smp(45,ix)+255000),findgroups(mcsec(ix)));...
        splitapply(@nanmean,min(189000,smp(65,ix)+255000),findgroups(mcsec(ix)))];
    out2 = out2/101972;
    
    
    
    %plotting
    xdk = figure;
    xv = 0.25:0.5:24;
    ss = [2,2,1,1];
    tt = {'','PHS','','SMS'};
    cc = [0,0,0.7,0.7];
    ll = {'-',':','-',':'};
    yy1 = [10e-9,10e-9,1.5e-10,1.5e-10];
    rr = {'(b)','','(a)',''};
    for i=1:4
        subplot(2,2,ss(i))
        plot(xv,out(i,:),'LineWidth',2,'LineStyle',ll{i},'Color',cc(i)*ones(1,3))
        xlim([0,24])
        ylim([0,yy1(i)])
        set(gca,'xtick',0:6:24)
        title(tt{i})
        hold on
        box off
        if i==1
            ylabel({'Layer 5 mean (modeled)';'Hydraulic Conductance (s-1)'})
        elseif i==3
            ylabel({'Layer 5 mean (implied)';'Hydraulic Conductance (s-1)'})
        end
        
        text(21,0.1*yy1(i),rr{i},'FontSize',14,'FontWeight','bold')
        
    end
    
    
    yy1 = [0.25,0.25,2,2];
    yy2 = [-0.15,-0.15,0,0];
    
    rr = {'(d)','','(c)',''};
    for i=1:4
        subplot(2,2,ss(i)+2)
        plot(xv,out2(i,:),'LineWidth',2,'LineStyle',ll{i},'Color',cc(i)*ones(1,3))
        xlim([0,24])
        set(gca,'xtick',0:6:24)
        xlabel('Hour of Day')
        ylim([yy2(i),yy1(i)])
        hold on
        box off
        ylabel({'Layer 5 mean';'\Delta\psi (MPa)'})
        if i==4
            legend('AMB','TFE','location','Southwest')
        end
        
        text(21,yy2(i)+0.92*(yy1(i)-yy2(i)),rr{i},'FontSize',14,'FontWeight','bold')
        
    end
    

    figure
    ix = year==2003;
    
    plot([0,4e-10],[0,4e-10],'k:')
    hold on
    plot(k(45,ix),k(65,ix),'.')
    xlim([0,4e-10])
    ylim([0,4e-10])
    
    xdk.Units = 'Inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    if hh(14)>1
        print(xdk,'../figs3/k','-dpdf')
    end
end

if hh(13)>0
    
    
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
    
    theq = zeros(4,4);
    cc = 0;
    for i = 1:2
        cc=cc+1;
        
        ix = vpd>1&vpd<1.0559&year>2001;
        max(fsds(ix))
        subplot('Position',[xx(i),0.53,0.4,0.41])
        qq = quantile(sssms(i,ix),[1/3,2/3]);
        qq = [-inf,qq,0];
        theq(cc,:) = qq;
        for j=[3,1,2]
            ix2 = ix&sssms(i,:)>=qq(j)&sssms(i,:)<qq(j+1);
            plot(fsds(ix2),btran(i+2,ix2),'.')
            xlim([0,1000])
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
        
        text(900,0.9,ss{i},'FontSize',14,'FontWeight','bold')
    end
    
    
    
    ssphs  = repmat(vegwp(4,mcsec==diurn(10)),48,1);
    ssphs  = ssphs(1:52551);
    ssphs2 = repmat(vegwp(8,mcsec==diurn(10)),48,1);
    ssphs2 = ssphs2(1:52551);
    ssphs  = [ssphs;ssphs2];
    
    
    
    ss = {'(c)','(d)'};
    for i = 1:2
        cc = cc+1;
        subplot('Position',[xx(i),0.08,0.4,0.41])
        ix = vpd>1&vpd<1.0559&year>2001;
        qq = quantile(ssphs(i,ix),[1/3,2/3]);
        qq = [-inf,qq,0];
        theq(cc,:) = qq;
        for j=[3,1,2]
            ix2 = ix&ssphs(i,:)>=qq(j)&ssphs(i,:)<qq(j+1);
            plot(fsds(ix2),btran(i,ix2),'.')
            xlim([0,1000])
            ylim([0,1])
            hold on
        end
        box off
        xlabel('Shortwave Radiation (W/m2)')
        if i==1
            legend('wettest','driest','intermediate','location','SouthWest')
            ylabel('Stress Factor (PHS)')
        else
            set(gca,'yticklabel',[])
            
        end
        text(900,0.9,ss{i},'FontSize',14,'FontWeight','bold')
        
    end
    
    
    xdk.Units = 'inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    if hh(13)>1
        print(xdk,'../figs3/suppstress','-dpdf')
    end
    
end


if hh(12)>0
    
    
    g = (year-2001)*365+doy;
    x1 = splitapply(@mean,zr*smp(41:60,:),g)';
    x2 = splitapply(@mean,zr*smp(61:80,:),g)';
    x3 = vegwp(4,mcsec==diurn(10))';
    x4 = vegwp(8,mcsec==diurn(10))';
    x  = [x1,x2,x3,x4]/101972;
    y1 = p(:,1);
    y2 = p(:,2);
    y3 = 1800*4e-7*splitapply(@sum,fctr',g');
    y  = [y1,y2,y3];
    
    xdk = figure;
    ix = p(:,1)>0;
    subplot(2,2,1)
    plot(x(ix,1),y(ix,5)-y(ix,1),'.')
    xlim([-4,0])
    ylim([-3,3])
    subplot(2,2,3)
    plot(x(ix,3),y(ix,3)-y(ix,1),'.')
    ylim([-3,3])
    xlim([-1,0])
    
    ix = p(:,2)>0;
    subplot(2,2,2)
    plot(x(ix,2),y(ix,6)-y(ix,2),'.')
    xlim([-4,0])
    ylim([-3,3])
    subplot(2,2,4)
    plot(x(ix,4),y(ix,4)-y(ix,2),'.')
    xlim([-1,0])
    ylim([-3,3])
    
    ss = {'SMSamb','SMStfe','PHSamb','PHStfe'};
    for i=1:4
        subplot(2,2,i)
        title(ss{i})
        ylabel({'Model-Obs';'Transpiration (mm/d)'})
        xlabel('Model Soil Potential (MPa)')
    end
    
    
    
    xdk.Units = 'Inches';
    xdk.Position = [2,2,7,5];
    xdk.PaperSize = [7,5];
    xdk.PaperPosition = [0,0,7,5];
    
    
    
    
    
end



if hh(11)>0
    g = (year-2001)*365+doy;
    x1 = splitapply(@mean,zr*smp(41:60,:),g)';
    x2 = splitapply(@mean,zr*smp(61:80,:),g)';
    x3 = vegwp(4,mcsec==diurn(10))';
    x4 = vegwp(8,mcsec==diurn(10))';
    x  = [x1,x2,x3,x4]/101972;
    y1 = p(:,1);
    y2 = p(:,2);
    y3 = 1800*4e-7*splitapply(@sum,fctr',g');
    y  = [y1,y2,y3];
    
    
    xdk = figure;
    c = 0;
    for j=[1,2]
        ix = p(:,j)>0;
        for i=[1,5]
            c = c+1;
            subplot(4,2,c)
            plot(x(ix,j),y(ix,i-1+j),'.','Color',[0.7,0.7,0.7])
            xlim([-4,0])
            ylim([0,6])
            
            if i==1
                ylabel({'Obs Sap Flux','(mm/d)'})
            else
                ylabel({'Model T','(mm/d)'})
            end
            set(gca,'ytick',0:3:6)
            
            if j==1
                if i==1
                    title('Observed Transpiration')
                else
                    title('Modeled Transpiration')
                end
            end
            
            
            
        end
    end
    
    xi = [3,3,4,4];
    yi = [1,3,2,4];
    
    for i=1:4
        ix = p(:,1+(i>2))>0;
        subplot(4,2,4+i)
        plot(x(ix,xi(i)),y(ix,yi(i)),'k.')
        xlim([-1,0])
        ylim([0,6])
        set(gca,'ytick',0:3:6)
        if i>2
            xlabel('Model Soil Potential (MPa)')
        end
        
        if i==1||i==3
            ylabel({'Obs Sap Flux','(mm/d)'})
        else
            ylabel({'Model T','(mm/d)'})
        end
        
    end
    
    ax1 = axes('Position',[0 0 1 1],'Visible','off');
    ss1 = {'SMS','SMS','PHS','PHS'};
    ss2 = {'amb','tfe','amb','tfe'};
    for i=1:4
        text(0.92,1.08-0.22*i,ss1{i},'FontSize',12)
        text(0.92,1.05-0.22*i,ss2{i},'FontSize',12)
    end
    
    xdk.Units = 'Inches';
    xdk.Position = [2,2,7,7];
    xdk.PaperSize = [7,7];
    xdk.PaperPosition = [0,0,7,7];
    
    
    if hh(11)>1
        print(xdk,'../figs3/transpiration_vs_smp','-dpdf')
    end
    
end







if hh(10)>0
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
    
    
    ssphs  = repmat(vegwp(4,mcsec==diurn(10)),48,1);
    ssphs  = ssphs(1:52551);
    ssphs2 = repmat(vegwp(8,mcsec==diurn(10)),48,1);
    ssphs2 = ssphs2(1:52551);
    ssphs  = [ssphs;ssphs2];
    
    
    n = sum(ix);
    plot(0:1/(n-1):1,sort(sssms(1,ix))/101972,'.')
    hold on
    plot(0:1/(n-1):1,sort(ssphs(1,ix))/101972,'.')
    plot([1/3,1/3],[-4,0],'k:')
    plot([2/3,2/3],[-4,0],'k:')
    ylim([-3,0])
    xlim([0,1])
    
    ix = ix&year>2001;
    n = sum(ix);
    plot(0:1/(n-1):1,sort(sssms(2,ix))/101972,'.')
    hold on
    plot(0:1/(n-1):1,sort(ssphs(2,ix))/101972,'.')
    plot([1/3,1/3],[-4,0],'k:')
    plot([2/3,2/3],[-4,0],'k:')
    ylim([-3,0])
    xlim([0,1])
    
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
    yy = [0.55,0.55,0.78,0.78];
    w  = 0.44;
    h  = 0.17;
    cc = [0,0,0.7,0.7];
    ss = '(c)(d)(a)(b)';
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
            l = legend('PHS');
            l.Position =[0.3135    0.5561    0.1210    0.0417];
            %l.Position =[0.32    0.7944    0.11    0.04];
        elseif i==4
            title('TFE')
        elseif i==3
            title('AMB')
            ll=legend('SMS');
            ll.Position =[0.3135    0.7861    0.1210    0.0417];
        end
        
        if i==2||i==4
            set(gca,'yticklabel',[])
        end
        text(10,0.2,ss((1:3)+(i-1)*3),'FontSize',14,'FontWeight','bold')
        box off
    end
    
    xf = 86400*1e-6*12; %umol/m2/s to g/m2/d
    yy = [0.09,0.09,0.32,0.32];
    ss = '(g)(h)(e)(f)';
    ix = year>2001;
    g = (year(ix)-2002)*365+doy(ix);
    out = xf*splitapply(@mean,fpsn(:,ix)',g');
    for i=1:4
        s=subplot('Position',[xx(i),yy(i),w,h]);
        plot([365,365],[0,10],'Color',[0.3,0.3,0.3])
        hold on
        
        plot(xv,out(:,i),'.','Color',cc(i)*ones(1,3))
        xlim([0,730])
        ylim([0,10])
        if i>2
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
    
    if hh(9)>1
        print(xdk,'../figs3/gpp','-dpdf')
    end
    
    disp('std of daily total gpp (gC/m2/d)')
    disp('PHSamb PHStfe SMSamb SMStfe')
    std(out)
    
    ix = year>2001;
    disp('total gpp kgC/m2')
    1800*sum(fpsn(:,ix),2)*1e-9*12
    
    ttest(fpsn(1,ix),fpsn(3,ix))
    ttest(fpsn(2,ix),fpsn(4,ix))
    
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
    
    theq = zeros(4,4);
    cc = 0;
    for i = 1:2
        cc=cc+1;
        
        ix = fsds>400&fsds<425&year>2001;
        subplot('Position',[xx(i),0.53,0.4,0.41])
        qq = quantile(sssms(i,ix),[1/3,2/3]);
        qq = [-inf,qq,0];
        theq(cc,:) = qq;
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
        xlim([0,3.2])
        set(gca,'xtick',0:3)
        set(gca,'xticklabel',[])
        text(2.87,0.1,ss{i},'FontSize',14,'FontWeight','bold')
    end
    
    
    
    ssphs  = repmat(vegwp(4,mcsec==diurn(10)),48,1);
    ssphs  = ssphs(1:52551);
    ssphs2 = repmat(vegwp(8,mcsec==diurn(10)),48,1);
    ssphs2 = ssphs2(1:52551);
    ssphs  = [ssphs;ssphs2];
    

            
    ss = {'(c)','(d)'};
    for i = 1:2
        cc = cc+1;
        subplot('Position',[xx(i),0.08,0.4,0.41])
        ix = fsds>400&fsds<425&year>2001;
        qq = quantile(ssphs(i,ix),[1/3,2/3]);
        qq = [-inf,qq,0];
        theq(cc,:) = qq;
        for j=[3,1,2]
            ix2 = ix&ssphs(i,:)>=qq(j)&ssphs(i,:)<qq(j+1);
            plot(vpd(ix2),btran(i,ix2),'.')
            ylim([0,1])
            hold on
        end
        box off
        xlabel('VPD (kPa)')
        if i==1
            legend('wettest','driest','intermediate','location','SouthWest')
            ylabel('Stress Factor (PHS)')
        else
            set(gca,'yticklabel',[])
            
        end
        text(2.87,0.1,ss{i},'FontSize',14,'FontWeight','bold')
        xlim([0,3.2])
        set(gca,'xtick',0:3)
        
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
    ss = '(e)(f)(c)(d)';
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
            text(0.2,5.4,ss((1:3)+(c-1)*3),'FontSize',14,'FontWeight','bold')
        end
    end
    
    ss = '(a)(b)';
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
        
        text(2/60*2+2002,5.4,ss((1:3)+(i-1)*3),'FontSize',14,'FontWeight','bold')
        
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
    
    
    ix = year==2003&month>8&month<12&mcsec>=diurn(25)&mcsec<=diurn(28);
    mean(btran(:,ix),2)
    
    
    if hh(2)>1
        print(xdk,'../figs3/fig4','-dpdf')
    end
    
    
end

if hh(56)>0
    g = (year-2001)*12+month;
    ix = mcsec>=diurn(25)&mcsec<=diurn(28);
    
    xv = cumsum(repmat(eomday(2001,1:12),1,3));
    xv = 2001+([0,xv(1:end-1)]+xv)/2/365;
    
    xdk = figure;
    ll = {'-',':'};
    
    for ee=[0,1]

        z = splitapply(@mean,vegwp(1+ee*4,ix),g(ix));
        set(gca,'ColorOrderIndex',1)
        plot(xv,z/101972,'LineWidth',2,'LineStyle',ll{ee+1})
        hold on

    end
    for ee=[0,1]
            y = splitapply(@mean,vegwp(4+ee*4,ix),g(ix));
                    set(gca,'ColorOrderIndex',2)
        plot(xv,y/101972,'LineWidth',2,'LineStyle',ll{ee+1})
    end
    ylim([-3,0])
    legend({'SunLeafAMB','SunLeafTFE','RootAMB','RootTFE'},'Location','Southwest') 
    
end

if hh(1)>0
    cols = [...
             0    0.4470    0.7410;...
    0.4940    0.1840    0.5560;...
    0.9290    0.6940    0.1250;...
    0.8500    0.3250    0.0980];
    %figure 2, SON2003 diurnal and monthly mean of veg water potential
    
    % calculate monthly mean midday lwp
    tt=25;
    ixd = mcsec>=diurn(tt)&mcsec<=diurn(tt+3);
    
    for i=1:2
        lwp(i,:) = 1/101972*splitapply(@mean,vegwp(1+(i-1)*4,ixd),month(ixd)+(year(ixd)-2001)*12);
        rwp(i,:) = 1/101972*splitapply(@mean,vegwp(4+(i-1)*4,ixd),month(ixd)+(year(ixd)-2001)*12);
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
    for i=4:-1:1
        plot(x,xf*out(i,:),'Color',cols(i,:),'LineWidth',2)
        hold on
    end
    xlim([0 24])
    ylim([-3 0])
    xlabel('Hour')
    ylabel({'Water Potential';'(MPa)'})
    title('AMB')
    text(1,-2.5,'(a)','FontSize',14,'FontWeight','bold')
    set(gca,'xtick',0:6:24)
    
    subplot('position',[0.545, 0.59, 0.425, 0.36])
    %plot(x,xf*out(5:8,:)')
    for i=4:-1:1
        plot(x,xf*out(i+4,:),'Color',cols(i,:),'LineWidth',2)
        hold on
    end
    set(gca,'xtick',0:6:24)
    ylim([-3 0])
    xlim([0 24])
    xlabel('Hour')
    title('TFE')
    set(gca,'yticklabel',[])
    l = legend('root','stem','shade-leaf','sun-leaf','location','southeast');
    l.Position(1) = 0.45;
    l.Position(2) = 0.65;
    set(gca,'xtick',0:6:24)
    text(21,-2.5,'(b)','FontSize',14,'FontWeight','bold')
    
    x = 2001+(0.5:36)/12;
    
    subplot('position',[0.1, 0.1, 0.87, 0.36])

    plot(x,lwp(1,:),'-','LineWidth',2)
    hold on
    set(gca,'ColorOrderIndex',1)
    plot(x,lwp(2,:),':','LineWidth',2)
    plot(x,rwp(1,:),'-','LineWidth',2)
    set(gca,'ColorOrderIndex',2)
    plot(x,rwp(2,:),':','LineWidth',2)
    ylim([-3,0])
    %xlim([0,36])
    set(gca,'ytick',-3:0.5:0)
    set(gca,'yticklabel',{'-3','','-2','','-1','',0})
    set(gca,'xtick',2001:0.5:2004)
    set(gca,'xticklabel',{'2001','','2002','','2003','','2004'})
    xlabel('Year')
    ylabel({'Midday Water';'Potential (MPa)'})
    
    %draw an arrow
    plot(2001+[10/12,10/12],[-2.6,-2.9],'k-','LineWidth',2)
    plot(2001+[10/12,9.7/12],[-2.88,-2.8],'k-','LineWidth',2)
    plot(2001+[10/12,10.3/12],[-2.88,-2.8],'k-','LineWidth',2)
    
    legend('AMB','TFE','Location','Southwest')
    x = (87-42.5)/42.5*24;
    text(2001+3*(21+x)/(24+x),-1-2.5/3*2,'(c)','FontSize',14,'FontWeight','bold')
    
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
        print(xdk,'../figs3/vwp','-dpdf')
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
    
    
    corr(out(ix1,1),amb(ix1))^2;
    corr(out(ix2,2),tfe(ix2))^2;
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

