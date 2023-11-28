c=1;  %c=1 production function; c=-1 cost function
T=5; %time periods
pmiss=10; % percentage of missing values
N=100; %spatial units
R=10; %number of replications 
NT=N*T;
beta=0.5;
rho=0.4;
delta=0.5;
mu=0.5;
sigmau=0.1;
sigmav=0.2;
par=[beta;rho;delta;mu;sigmau;sigmav];

%import spatial data
spatialinfo=readtable("N100.csv"); %to be changed depending on N
lonlat=table2array(spatialinfo);
lon=lonlat(:,1);
lat=lonlat(:,2);
%create contiguity matrix
M=xy2cont(lon,lat);
M2=full(M);
W=sparse(kron(eye(T),M2)); 

%find eigenvalues min,max
%opt.tol=1e-3; opt.disp=0;
%lambda=eigs(sparse(W),speye(NT),1,'SR',opt);  
%rmin=1/lambda;   
%rmax=1;

%use Pace and Barry, 1999 MC approximation
%order=50; iter=30; 
%out=lndetmc(order,iter,W,rmin,rmax);
%tt=rmin:.001:rmax; 
%outi=interp1(out.rho,out.lndet,tt','spline');
%detval=[tt' outi];

parhat=zeros(6,R);

for i=1:N
    fe(:,i)=unifrnd(0,1);
    x(:,i)=normrnd(fe(:,i),1,T,1);
    z(:,i)=normrnd(0,1,T,1);
end

h=exp(z*delta);

for j=1:R
    
    for i=1:N
    v(:,i)=normrnd(0,sqrt(sigmav),T,1);
    distribu=makedist('Normal','mu',exp(mu),'sigma',sqrt(sigmau));
    distributrunc=truncate(distribu,0,inf);
    u=random(distributrunc,N,1)';
    end
    
    X=reshape(x',NT,1);
    alfa=repmat(fe',T,1);
    V=reshape(v',NT,1);
    U=reshape(h'.*u',NT,1);
    t=(eye(NT)-rho*W);
    Y=t\alfa+t\(X*beta)+t\V-c*t\U;
    
    %set data as missing at random
    x2=x;
    for i=1:T
    x2(i,randperm(N,pmiss/100*N))=missing;
    m=ismissing(x2);
    end
    
    z2=z;
    z2(m)=missing;

    dx=x2-repmat(mean(x2,'omitNaN'),T,1);
    dz=z2-repmat(mean(z2,'omitNaN'),T,1);
    dh=exp(dz*delta);
    dv=v-repmat(mean(v),T,1);
    du=dh.*u;

    y=reshape(Y',N,T)';
    y(m)=missing;
    Y2=reshape(y',NT,1);
    Y22=Y2;
    Y22(isnan(Y2))=0;
    neigh=W*Y22;
    neigh(isnan(Y2))=NaN;
    neigh=reshape(neigh,N,T)';
    dneigh=neigh-repmat(mean(neigh,'omitNaN'),T,1);
    dy=rho*dneigh+beta*dx+dv-c*du; 

    beta0=0.2;
    rho0=0.6;
    delta0=0.7;
    mu0=0.4;
    sigmau0=0.1;
    sigmav0=0.3;
    par0=[beta0;rho0;delta0;mu0;sigmau0;sigmav0];
    fun=@(parhat)loglik(parhat,dy,dx,dz,N,T,dneigh,c,W);
    
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    lb=[-inf,-1,-inf,-inf,0,0];
    ub=[inf,1,inf,inf,inf,inf];
    nonlncon=[];
    options=optimoptions(@fmincon,'MaxFunEval',1e17,'MaxIter',1e17,'OptimalityTolerance',1e-20, 'StepTol', 1e-100);
    [parhat(:,j)]=fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlncon,options);

    j

end

diff=parhat-par; 
stime=mean(parhat,2);
bias=mean(diff,2);
diff2=(diff).^2;
mse=mean(diff2,2);

par2=[stime(1,:);stime(2,:);stime(3,:);stime(4,:);stime(5,:);stime(6,:)];
options2=optimoptions(@fminunc,'MaxFunEval',1e17,'MaxIter',1e17,'OptimalityTolerance',1e-17,'StepTol',1e-100);
[p1,fval,exitflag,output,grad,hessian1]=fminunc(fun,par2,options2);
SD=sqrt(diag(inv(hessian1)));

