
clear
close all

%ignore this part
%I just used it to generate random numbers for p and n
%However, your variables p and n should look similar.
%They both need to be in order and the same length
n = 10:10:300;
pp = 0.97;
for i=1:30
p(i) = pp;
pp = rand*pp/5+4*pp/5;
end


%this part I plot up the points
plot(n,p,'rx')


%you can estimate where the two missing points might be
p1 = 1;
n1 = (1-p(1))/p(1)*n(1);
p2 = 0;
n2 = n(end)+ p(end)/(p(end-1)-p(end))*(n(end)-n(end-1));

%then i can plot up the curve
hold on
plot([n1,n,n2],[p1,p,p2],'b-')
xlim([0,1.05*n2])

%here is where I compute the area


main_area = sum(0.5*(p(2:end)+p(1:end-1)).*(n(2:end)-n(1:end-1)))

%first_bit = ???
%last_bit = ???
