function [ a ] = getmore( filelist,offset,n,varlist,vard )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ne  = length(filelist);

for vv=1:length(varlist)
    out=zeros(sum(vard(vv,:)),n);
    for ee=1:ne
        if vard(vv,ee)>0
            temp  = ncread(filelist{ee},varlist{vv});
            if ee>1
                idx   = sum(vard(vv,1:ee-1))+(1:vard(vv,ee));
            else
                idx   = (1:vard(vv,ee));
            end
            if vard(vv,ee)>1
                
                out(idx,:)=temp(1,1:vard(vv,ee),1+offset:end);
            else
                out(idx,:)=temp(1,1+offset:end);
            end
        end
    end
    
    assignin('caller',lower(varlist{vv}),out)
end



a=1;
end

