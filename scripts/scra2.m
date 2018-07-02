clear
close all

ix  = zeros(40,80);
ix(:) = 1:3200;



out = zeros(40,80);
for i=1:40
    for j=1:80
        out(i,j) = sqrt(i)*sqrt(j);
    end
end

p = ix(out>33);
r = 1+floor(length(p)*rand);

j = 1+floor(p(r)/40.001);
i = mod(p(r),40);
if i==0
    i=40;
end


xdk = figure;

%out2 = min(out,2000);
out2 = min(out,33);
out2(isnan(out))= nan;
out2(i,j) = out(i,j)*2;

subplot(1,2,1)
imagesc(out',[0,55])
set(gca,'YDir','Normal')
c = colorbar;

subplot(1,2,2)
imagesc(out2',[0,55])
set(gca,'YDir','Normal')
c = colorbar;

xdk.Position=[360   472   858   226];