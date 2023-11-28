T=5;
N=100;
R=5; 
NT=N*T;
beta=0.5;
rho=0.3;
theta=0.3;
phi=0.2;
delta=0.5;
sigmau=0.1;
sigmav=0.1;
sigma2=sigmau+sigmav;
lambda=sigmau/sigma2;
par=[beta;rho;theta;phi;delta;sigma2;lambda];

%Spatial weight matrix W
N100=readtable('N100');
lonlat=table2array(N100);
lon=lonlat(:,1);
lat=lonlat(:,2);
%Binary Contiguity 
M=xy2cont(lon,lat);
M2=full(M);
W=kron(eye(T),M2);
W=sparse(W);
%Inverse Distance
%M1=pdweight(lon,lat,0,15000,1);
%M3=full(M1);
%W1=kron(eye(T-1),M3);
%W1=sparse(W1);


%Pace and Barry 1999 MC approximation to compute log-determinant 
opt.tol=1e-3; opt.disp=0;
eig=eigs(M,speye(N),1,'SR',opt);  
rmin=1/eig;   
rmax=1;
order=50; 
iter=30; 
out=lndetmc(order,iter,M,rmin,rmax);
tt=rmin:.001:rmax; 
outi=interp1(out.rho,out.lndet,tt','spline');
detval=[tt' outi];

X=normrnd(0,1,NT,1);
Z=normrnd(0,1,NT,1);
parhat=zeros(7,R);
I=eye(NT);
t=(I-rho*W);

for j=1:R
   
    v=(sqrt(sigmav))*randn(NT,1);
    distribw=makedist('Normal','mu',0,'sigma',sqrt(sigmau));
    w=zeros(NT,1);
    for i=1:NT
    mu=Z*phi+W*Z*delta;
    distribwtrunc(i)=truncate(distribw,-mu(i),inf);    
    w(i)=random(distribwtrunc(i));    
    w=transpose(w);
    end
    
    y=t\(X*beta)+t\(W*X*theta)+t\v-t\(Z*phi+W*Z*delta+w); 
    
    beta0=0.87;
    rho0=0.10;
    theta0=0.40;
    phi0=0.44;
    delta0=0.61;
    sigma20=0.33;
    lambda0=0.34;
    par0=[beta0;rho0;theta0;phi0;delta0;sigma20;lambda0];
    fun=@(parhat)loglik(parhat,y,X,Z,W,M,NT,N,T);
    
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    lb=[-inf,-1,-inf,-inf,-inf,1e-15,1e-15];
    ub=[inf,1,inf,inf,inf,inf,0.99];
    nonlncon=[];
    options = optimoptions(@fmincon,'MaxFunEval',1e17,'MaxIter',1e17,'OptimalityTolerance',1e-17, 'StepTol', 1e-100);
    [parhat(:,j)]=fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlncon,options);
   
end 

est=mean(parhat,2);
diff=parhat-par; 
bias=mean(diff,2);
diff2=(diff).^2;
MSE=mean(diff2,2);

%Wang and Ho (2010) function to compute robust variance covariance...
    %matrix numerically
[VCVROBUST,~,~,~,hess]=robustvcv(@loglik,est,0,y,X,Z,W,M,NT,N,T);
VCV=hess^(-1)/(NT);
ste=sqrt(diag(VCVROBUST));
robste=sqrt(diag(VCV));
tstat=est./robste;
pvalue=tpdf(tstat,N*(T-1));

