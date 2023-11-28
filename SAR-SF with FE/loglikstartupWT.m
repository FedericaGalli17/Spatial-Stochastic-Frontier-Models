function [val,l]=loglikstartupWT(parhat,dY,dL,dK,N,dneigh,c,detval,dbrev,dlev,dint,dsize,dliquid,dred,dLL,dKK,dLK)

gsize=detval(2,1)-detval(1,1);
i1=find(detval(:,1)<=parhat(6)+gsize);
i2=find(detval(:,1)<=parhat(6)-gsize);
i1=max(i1);
i2=max(i2);
index=round((i1+i2)/2);
if isempty(index)
index=1;
end

detm=detval(index,2); 

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
muast=(parhat(13)/parhat(14)-c*transpose(deps)*pinv(sigma)*dh)/(transpose(dh)*pinv(sigma)*dh+1/parhat(14));
sast=1/(transpose(dh)*pinv(sigma)*dh+1/parhat(14));

l(:,j)=-0.5*(t-1)*log(2*pi)-0.5*(t-1)*log(parhat(15))-0.5*transpose(deps)*pinv(sigma)*deps+0.5*((muast^2/sast)-(parhat(13)^2/parhat(14)))+log(sqrt(sast)*normcdf(muast/sqrt(sast)))-log(sqrt(parhat(14))*normcdf(parhat(13)/sqrt(parhat(14))));              

end


val=-detm-sum(l);

end



