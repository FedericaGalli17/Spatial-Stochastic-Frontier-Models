function [sast,mu]=scores(parhat,dY,dL,dK,N,dneigh,c,dbrev,dlev,dint,dsize,dliquid,dred,dLL,dKK,dLK)

for j=1:N

dz1=rmmissing(dsize(:,j));
dz2=rmmissing(dbrev(:,j));
dz3=rmmissing(dlev(:,j));
dz4=rmmissing(dliquid(:,j));
dz5=rmmissing(dint(:,j));
dz6=rmmissing(dred(:,j));
dy=rmmissing(dY(:,j));
dl=rmmissing(dL(:,j));
dk=rmmissing(dK(:,j));
dn=rmmissing(dneigh(:,j));
dll=rmmissing(dLL(:,j));
dkk=rmmissing(dKK(:,j));
dlk=rmmissing(dLK(:,j));

t=sum(~isnan(dY(:,j)));
    
d=(1-1/t)*parhat(15)*ones(t,1);
off=-1/t*parhat(15)*ones(t);
off1=triu(off,1);
sigma=off1+diag(d)+off1';    

dh=exp(dz1*parhat(7)+dz2*parhat(8)+dz3*parhat(9)+dz4*parhat(10)+dz5*parhat(11)+dz6*parhat(12));
deps=dy-parhat(1)*dl-parhat(2)*dk-parhat(3)*dll-parhat(4)*dkk-parhat(5)*dlk-parhat(6)*dn;

mu(:,j)=(parhat(13)/parhat(14)-c*transpose(deps)*pinv(sigma)*dh)/(transpose(dh)*pinv(sigma)*dh+1/parhat(14));
sast(:,j)=1/(transpose(dh)*pinv(sigma)*dh+1/parhat(14));

end 


end 

