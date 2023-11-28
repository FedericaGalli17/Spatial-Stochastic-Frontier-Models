c=1; %c=1 for production function and c=-1 for cost function
T=5;

dati=readtable('datiSU.csv'); %import data
M=readtable('mat50.txt'); %import spatial weight matrix
dati=table2array(dati);
M=table2array(M);
M(:,1)=[];
N=length(dati);

Y=dati(:,1:5);
L=dati(:,6:10);
LL=dati(:,6:10).^2;
K=dati(:,11:15);
KK=dati(:,11:15).^2;
LK=L.*K;
size=dati(:,16:20);
brev=dati(:,21:25);
lev=dati(:,26:30);
liquid=dati(:,31:35);
int=dati(:,36:40);
red=dati(:,41:45);

%drop missing data in spatial weight matrices
%t=1
M1=M;
indices=isnan(Y(:,1));
M1(:,indices)=0; 
M1(indices,:)=0;  
W1=M;
W1(:,indices)=[]; 
W1(indices,:)=[];  
%t=2
M2=M;
indices=isnan(Y(:,2));
M2(indices,:)=0;  
M2(:,indices)=0; 
W2=M;
W2(:,indices)=[]; 
W2(indices,:)=[];  
%t=3
M3=M;
indices=isnan(Y(:,3));
M3(indices,:)=0; 
M3(:,indices)=0; 
W3=M;
W3(:,indices)=[]; 
W3(indices,:)=[];  
%t=4
M4=M;
indices=isnan(Y(:,4));
M4(indices,:)=0; 
M4(:,indices)=0; 
W4=M;
W4(:,indices)=[]; 
W4(indices,:)=[];  
%t=5
M5=M;
indices=isnan(Y(:,5));
M5(indices,:)=0; 
M5(:,indices)=0; 
W5=M;
W5(:,indices)=[]; 
W5(indices,:)=[];  

%row-normalize matrices
M1=M1./sum(M1,2);
M1(isnan(M1))=0;
M2=M2./sum(M2,2);
M2(isnan(M2))=0;
M3=M3./sum(M3,2);
M3(isnan(M3))=0;
M4=M4./sum(M4,2);
M4(isnan(M4))=0;
M5=M5./sum(M5,2);
M5(isnan(M5))=0;
W1=W1./sum(W1,2);
W1(isnan(W1))=0;
W2=W2./sum(W2,2);
W2(isnan(W2))=0;
W3=W3./sum(W3,2);
W3(isnan(W3))=0;
W4=W4./sum(W4,2);
W4(isnan(W4))=0;
W5=W5./sum(W5,2);
W5(isnan(W5))=0;

%block-diagonal matrices
W=blkdiag(M1,M2,M3,M4,M5);
Wnan=blkdiag(W1,W2,W3,W4,W5);

%within-transposed data
dY=transpose(Y)-repmat(mean(transpose(Y),'omitNaN'),T,1);
dL=transpose(L)-repmat(mean(transpose(L),'omitNaN'),T,1);
dLL=transpose(LL)-repmat(mean(transpose(LL),'omitNaN'),T,1);
dK=transpose(K)-repmat(mean(transpose(K),'omitNaN'),T,1);
dKK=transpose(KK)-repmat(mean(transpose(KK),'omitNaN'),T,1);
dLK=transpose(LK)-repmat(mean(transpose(LK),'omitNaN'),T,1);
dsize=transpose(size)-repmat(mean(transpose(size),'omitNaN'),T,1);
dbrev=transpose(brev)-repmat(mean(transpose(brev),'omitNaN'),T,1);
dlev=transpose(lev)-repmat(mean(transpose(lev),'omitNaN'),T,1);
dliquid=transpose(liquid)-repmat(mean(transpose(liquid),'omitNaN'),T,1);
dint=transpose(int)-repmat(mean(transpose(int),'omitNaN'),T,1);
dred=transpose(red)-repmat(mean(transpose(red),'omitNaN'),T,1);

%spatial lag Y   
Y2=vertcat(Y(:,1),Y(:,2),Y(:,3),Y(:,4),Y(:,5));
Y22=Y2;
Y22(isnan(Y2))=0;
neigh=W*Y22;
neigh(isnan(Y2))=NaN;
neigh=[neigh(1:N,:) neigh(N+1:2*N,:) neigh(2*N+1:3*N,:) neigh(3*N+1:4*N,:) neigh(4*N+1:5*N,:)];
dneigh=transpose(neigh)-repmat(mean(transpose(neigh),'omitNaN'),T,1);

%find eigenvalues min,max
opt.tol=1e-3; opt.disp=0;
lambda=eigs(sparse(W),speye(N*T),1,'SR',opt);  
rmin=1/lambda;   
rmax=1;

%use Pace and Barry, 1999 MC approximation
order=50; iter=30; 
out=lndetmc(order,iter,W,rmin,rmax);
tt=rmin:.001:rmax; 
outi=interp1(out.rho,out.lndet,tt','spline');
detval=[tt' outi];

par0=[0.5,0.3,0.3,0.1,0.1,0.1,0.2,0.1,-0.1,-0.1,-0.1,-0.1,0,0.5,0.5];
fun=@(parhat)loglikstartupWT(parhat,dY,dL,dK,N,dneigh,c,detval,dbrev,dlev,dint,dsize,dliquid,dred,dLL,dKK,dLK);    
A=[];
b=[];
Aeq=[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];
beq=[0];
lb=[-inf,-inf,-inf,-inf,-inf,-1,-inf,-inf,-inf,-inf,-inf,-inf,-inf,0,0];
ub=[inf,inf,inf,inf,inf,1,inf,inf,inf,inf,inf,inf,inf,inf,inf];
nonlncon=[];
options=optimoptions(@fmincon,'Display','iter','MaxFunEval',1e17,'MaxIter',1e17,'OptimalityTolerance',1e-17, 'StepTol', 1e-100);
[parhat,val]=fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlncon,options);
parhat %parameter estimates

%SD
[VCVROBUST,~,~,~,hess]=robustvcv(@loglikstartupWT,parhat,0,dY,dL,dK,N,dneigh,c,detval,dbrev,dlev,dint,dsize,dliquid,dred,dLL,dKK,dLK);
VCV=hess^(-1)/(N);
ste=sqrt(diag(VCVROBUST));
robste=sqrt(diag(VCV));
tstat=parhat./robste';
pvalue=tpdf(tstat,N*T-1);

%efficiency scores 
[sast,mu]=scores(parhat,dY,dL,dK,N,dneigh,c,dbrev,dlev,dint,dsize,dliquid,dred,dLL,dKK,dLK);
hi=exp(size*parhat(7)+brev*parhat(8)+lev*parhat(9)+liquid*parhat(10)+int*parhat(11)+red*parhat(12));
h=hi';
va=h.*(mu+(normpdf(mu./sqrt(sast)).*sqrt(sast))./(normcdf(mu./sqrt(sast))));
TE1=exp(-va);
TE_mean=mean(TE1,'omitNaN');

%marginal effects
elL=parhat(1)+parhat(3)*L+parhat(5)*K;
eL=mean(elL,'omitNaN');
eL=mean(eL);
elK=parhat(2)+parhat(4)*K+parhat(5)*L;
eK=mean(elK,'omitNaN');
eK=mean(eK);

elL1=rmmissing(elL(:,1));
elL2=rmmissing(elL(:,2));
elL3=rmmissing(elL(:,3));
elL4=rmmissing(elL(:,4));
elL5=rmmissing(elL(:,5));
elaL=vertcat(elL1,elL2,elL3,elL4,elL5);
elK1=rmmissing(elK(:,1));
elK2=rmmissing(elK(:,2));
elK3=rmmissing(elK(:,3));
elK4=rmmissing(elK(:,4));
elK5=rmmissing(elK(:,5));
elaK=vertcat(elK1,elK2,elK3,elK4,elK5);

I=eye(length(Wnan));
t=(I-parhat(6)*Wnan)^(-1);

dY_dL=t*diag(elaL);
dY_dK=t*diag(elaK);

Ldir=mean(diag(dY_dL)); %direct effect L
Ltot=mean(sum(dY_dL)); %total effect L
Lind=Ltot-Ldir; %indirect effect L
Kdir=mean(diag(dY_dK)); %direct effect K
Ktot=mean(sum(dY_dK)); %total effect K  
Kind=Ktot-Kdir; %indirect effect K

%delta method
lav=[L(:,1);L(:,2);L(:,3);L(:,4);L(:,5)];
lav=rmmissing(lav);
kap=[K(:,1);K(:,2);K(:,3);K(:,4);K(:,5)];
kap=rmmissing(kap);

%L
d1=(I-parhat(:,6)*Wnan)^(-1);
d2=(I-parhat(:,6)*Wnan)^(-1)*(diag(lav));
d3=(I-parhat(:,6)*Wnan)^(-1)*(diag(kap));
d4=(I-parhat(:,6)*Wnan)^(-2)*Wnan*(diag(elaL));
dir1=mean(diag(d1));
dir2=mean(diag(d2));
dir3=mean(diag(d3));
dir4=mean(diag(d4));
tot1=mean(sum(d1));
tot2=mean(sum(d2));
tot3=mean(sum(d3));
tot4=mean(sum(d4));
ind1=tot1-dir1;
ind2=tot2-dir2;
ind3=tot3-dir3;
ind4=tot4-dir4;
d=[dir1,dir2,dir3,dir4]';
vmat=VCV([1 3 5 6],[1 3 5 6]);
sd=d'*vmat*d;
sddirL=sqrt(sd); %sd direct effect
d=[ind1,ind2,ind3,ind4]';
vmat=VCV([1 3 5 6],[1 3 5 6]);
sd=d'*vmat*d;
sdindL=sqrt(sd); %sd indirect effect
d=[tot1,tot2,tot3,tot4]';
vmat=VCV([1 3 5 6],[1 3 5 6]);
sd=d'*vmat*d;
sdtotL=sqrt(sd); %sd total effect
td=Ldir/sddirL;
ti=Lind/sdindL;
tt=Ltot/sdtotL;
pvalued=tpdf(td,N*T-1);
pvaluei=tpdf(ti,N*T-1);
pvaluet=tpdf(tt,N*T-1);

%K
d1=(I-parhat(:,6)*Wnan)^(-1);
d2=(I-parhat(:,6)*Wnan)^(-1)*(diag(kap));
d3=(I-parhat(:,6)*Wnan)^(-1)*(diag(lav));
d4=(I-parhat(:,6)*Wnan)^(-2)*Wnan*(diag(elaK));
dir1=mean(diag(d1));
dir2=mean(diag(d2));
dir3=mean(diag(d3));
dir4=mean(diag(d4));
tot1=mean(sum(d1));
tot2=mean(sum(d2));
tot3=mean(sum(d3));
tot4=mean(sum(d4));
ind1=tot1-dir1;
ind2=tot2-dir2;
ind3=tot3-dir3;
ind4=tot4-dir4;
d=[dir1,dir2,dir3,dir4]';
vmat=VCV([2 4 5 6],[2 4 5 6]);
sd=d'*vmat*d;
sddirK=sqrt(sd); %sd direct effect
d=[ind1,ind2,ind3,ind4]';
vmat=VCV([2 4 5 6],[2 4 5 6]);
sd=d'*vmat*d;
sdindK=sqrt(sd); %sd indirect effect
d=[tot1,tot2,tot3,tot4]';
vmat=VCV([2 4 5 6],[2 4 5 6]);
sd=d'*vmat*d;
sdtotK=sqrt(sd); %sd total effect
td=Kdir/sddirK;
ti=Kind/sdindK;
tt=Ktot/sdtotK;
pvalued=tpdf(td,N*T-1);
pvaluei=tpdf(ti,N*T-1);
pvaluet=tpdf(tt,N*T-1);



