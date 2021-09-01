function [ f ] = calfactor(x,m,gbx) 
d=zeros(1,m); 
dg=0; 
for i=1:m 
    for j=1:m 
        d(i)=d(i)+pdist([x(i,:);x(j,:)]); 
    end 
    d(i)=d(i)/(m-1); 
end 
for j=1:m 
    dg=dg+pdist([gbx;x(j,:)]); 
end 
dg=dg/(m-1); 
d=[d dg]; 
f=(dg-min(d))/(max(d)-min(d)); 
end 