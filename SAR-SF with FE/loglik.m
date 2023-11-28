function [val]=loglik(parhat,dy,dx,dz,N,T,dneigh,c,W)

 if parhat(5)<=0
    parhat(5)=0.5;
 end    
 if parhat(6)<=0
    parhat(6)=0.5;
 end

idx=ismissing(reshape(dy',N*T,1));
W(idx,:)=[];
W(:,idx)=[];
nm=sum(idx);
detm=log(det(eye(N*T-nm)-parhat(2)*W)); 


for j=1:N 

dZ=rmmissing(dz(:,j));
dY=rmmissing(dy(:,j));
dX=rmmissing(dx(:,j));
dn=rmmissing(dneigh(:,j));
t=sum(~isnan(dy(:,j)));    

d=(1-1/t)*parhat(6)*ones(t,1);
off=-1/t*parhat(6)*ones(t);
off1=triu(off,1);
sigma=off1+diag(d)+off1';    

dh=exp(dZ*parhat(3));
deps=dY-parhat(2)*dn-parhat(1)*dX;
muast=(exp(parhat(4))./parhat(5)-c*transpose(deps)*pinv(sigma)*dh)/(transpose(dh)*pinv(sigma)*dh+1./parhat(5));
sast=1/(transpose(dh)*pinv(sigma)*dh+1./parhat(5));

l(:,j)=-0.5*(t-1)*log(2*pi)-0.5*(t-1)*log(parhat(6))-0.5*transpose(deps)*pinv(sigma)*deps+0.5*((muast^2/sast)-(exp(parhat(4))^2/parhat(5)))+log(sqrt(sast)*normcdf(muast/sqrt(sast)))-log(sqrt(parhat(5))*normcdf(exp(parhat(4))/sqrt(parhat(5))));              


end

val=-detm-sum(l);

end


