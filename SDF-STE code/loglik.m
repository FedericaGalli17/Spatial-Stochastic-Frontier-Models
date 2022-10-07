function [val,ll]= loglik(parhat,y,X,Z,W,M,NT,N,T)
    
    %For speed use the Pace and Barry (1999) approximation instead of the
        %direct computation of the log-determinant 
    %gsize=detval(2,1)-detval(1,1);
    %i1=find(detval(:,1)<=parhat(2)+gsize);
    %i2=find(detval(:,1)<=parhat(2)-gsize);
    %i1=max(i1);
    %i2=max(i2);
    %index=round((i1+i2)/2);
    %if isempty(index)
    %index=1;
    %end
    %detm=detval(index,2); 

    I=eye(N);
    detm=log(det(I-parhat(2)*M));
    
    var1=parhat(6);
    var2=parhat(7);
    den1=sqrt(var1*var2);
    eff=Z*parhat(4)+W*Z*parhat(5); 
    mod=y-X*parhat(1)-W*y*parhat(2)-W*X*parhat(3); 
    den2=sqrt(var1*var2*(1-var2));
    k=(eff*(1-var2))-(mod*var2); 
    
    val=T*detm-(NT/2)*(log(var1)+log(2*pi))-(1/2)*sum(((eff+mod).^2)./var1)-sum(log(normcdf(eff./den1))-log(normcdf(k./den2)));
    val=-val;
    
    ll=detm-(T/2)*(log(var1)+log(2*pi))-(1/2)*(((eff+mod).^2)./var1)-(log(normcdf(eff./den1))-log(normcdf(k./den2)));
    
end













