T=5;
N=300;
R=1000; 
NT=N*T;
beta=0.5;
rho=0.3;
theta=0.3;
delta=0.5;
tau=0.3;
phi=0.3;
sigma2u=0.1;
sigma2v=0.2;
par=[beta;rho;theta;delta;tau;phi;sigma2u;sigma2v];

%Import Geographic Coordinates
coord=readtable("lonlat");
lonlat=table2array(coord);
lon=lonlat(:,1);
lat=lonlat(:,2);
%Binary Contiguity W
M=xy2cont(lon,lat);
W=sparse(M);
%Inverse Distance W
%Change the upper bound distance (1500) to set a truncation point
M1=pdweight(lon,lat,0,1500,1);
W1=sparse(M1);

%Find min-max eigenvalues
opt.tol=1e-3; opt.disp=0;
lambda=eigs(sparse(W),speye(N),1,'SR',opt);  
rmin=1/lambda;   
rmax=1;

%Use Pace and Barry (1999) MC approximation
order=50; iter=30; 
out=lndetmc(order,iter,W,rmin,rmax);
tt=rmin:.001:rmax; 
outi=interp1(out.rho,out.lndet,tt','spline');
detval=[tt' outi];

%Generate Data
for k=1:T
    varname=['X',num2str(k)];
    assignin('caller',varname,normrnd(0,1,N,1))
end

for k=1:T
    varname=['Z',num2str(k)];
    assignin('caller',varname,normrnd(0,1,N,1))
end

H1=exp(Z1*delta);
H2=exp(Z2*delta);
H3=exp(Z3*delta);
H4=exp(Z4*delta);
H5=exp(Z5*delta);

parhat=zeros(8,R);

for j=1:R
    K=(eye(N)-phi*M)^(-1);
    matvar=sigma2v*K*transpose(K); 
    mu=zeros(N,1);
    v1=mvnrnd(mu,matvar)';
    v2=mvnrnd(mu,matvar)';
    v3=mvnrnd(mu,matvar)';
    v4=mvnrnd(mu,matvar)';
    v5=mvnrnd(mu,matvar)';

    distribu=makedist('Normal','mu',0,'sigma',sqrt(sigma2u));
    distributrunc=truncate(distribu,0,inf);
    u1=random(distributrunc);
    u2=random(distributrunc);
    u3=random(distributrunc);
    u4=random(distributrunc);
    u5=random(distributrunc);
    	    
    I=eye(N);
    t=(I-rho*M);
    y1=t\(X1*beta)+t\(M*X1*theta)+t\v1-t\((I-tau*M)\H1*u1); 
    y2=t\(X2*beta)+t\(M*X2*theta)+t\v2-t\((I-tau*M)\H2*u2); 
    y3=t\(X3*beta)+t\(M*X3*theta)+t\v3-t\((I-tau*M)\H3*u3); 
    y4=t\(X4*beta)+t\(M*X4*theta)+t\v4-t\((I-tau*M)\H4*u4); 
    y5=t\(X5*beta)+t\(M*X5*theta)+t\v5-t\((I-tau*M)\H5*u5); 
    
    beta0=0.5;
    rho0=0.5;
    theta0=0.5;
    delta0=0.5;
    tau0=0.5;
    phi0=0.5;
    sigmau20=0.7;
    sigmav20=0.4;
    par0=[beta0;rho0;theta0;delta0;tau0;phi0;sigmau20;sigmav20];
    fun=@(parhat)loglik2_T5(parhat,y1,y2,y3,y4,y5,X1,X2,X3,X4,X5,...    
        Z1,Z2,Z3,Z4,Z5,M,N,detval);
    
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    lb=[-inf,-1,-inf,-inf,-1,-1,0,0];
    ub=[inf,1,inf,inf,1,1,inf,inf];
    nonlncon=[];
    options=optimoptions(@fmincon,'MaxFunEval',1e17,'MaxIter',1e17,...
        'OptimalityTolerance',1e-17,'StepTol',1e-17);
    [parhat(:,j)]=fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlncon,options);
   
end

diff=parhat-par; 
est=mean(parhat,2);
bias=mean(diff,2);
diff2=(diff).^2;
MSE=mean(diff2,2);

%Wang and Ho (2010) function to compute robust variance covariance...
    %matrix numerically
[VCVROBUST,~,~,~,hess]=robustvcv(@loglik2_T5,est,0,y1,y2,y3,y4,y5,...
X1,X2,X3,X4,X5,Z1,Z2,Z3,Z4,Z5,W,N,detval);
VCV=hess^(-1)/(NT);
ste=sqrt(diag(VCVROBUST));
robste=sqrt(diag(VCV));
tstat=est./robste;
pvalue=tpdf(tstat,N*(T-1));

