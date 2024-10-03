function [output] = PBVARX_HIERARCHICAL_FUN_COMP(Yraw,Xraw,Iraw,varparaminputs)

% This function estimates the partially pooled VAR  

randn('seed',2); %#ok<RAND>
rand('seed',2);  %#ok<RAND>

% General notes: ALPHA refers to the matrix of all coefficients and alpha
% to the vectorised form. BETA refers to a matrix of VAR slope
% coefficients, beta to the vector of coefficients, GAMMA to the matrix of
% coefficients on the deterministic components and gamma to the
% corresponding vector. SIGMA is the covariance matrix of the residuals.
% BETAbar is the common cross-country mean of slopes, 

% In terms of hyperparameters, l refers to the scaling covariance matrix of
% the distribution of slope parameters, while v and s are the degrees of
% freedom and the shape parameters in the prior iG distribution of l. L
% refers to the unscaled covariance matrix of beta.

% A "c" at the end of any of these variables implies it is country
% specific.


% Written by Saleem Bahaj (Contact: saleembahaj@googlemail.com)
% Last updated: 28th December 2018

% Note: 
% This function is built upon BVAR code by Koop and Korobolis and relies
% upon auxiliary functions coded by Canova and Gambetti.


%% ========================== PRELIMINARIES ===============================
% ------------------------------ Priors -----------------------------------

% Prior for the reduced form BVAR model is a modified version of Jarocinski (2010).
% 1. betac                ~ Normal(betabar,l1*L1c)    - drawn from common cross country distibution (shrinkage) 
% 2. p(betabar),p(gammac) ~ 1                         - non-informative 
% 3. p(SIGMAc)            ~ iW(S,vbar)                - drawn from common cross country distibution (shrinkage)
% 4. p(S)                 ~ det(S)^-0.5(M+1)          - diffuse
% 5. l1                    ~ iG(v,s)                  - see Gelman et al (2003), Gelman (2006)
% 6. L1c is non-stochastic and calculated from the data.
% 7. Cross country covariance matrix is S/vbar
% Prior for the identification model is along a similar lines:
% 1. p(upsilonc)          ~ Normal(upsilonbar,l2*L2c) - drawn from common cross country distibution (shrinkage) 
% 2. p(upsilonbar)        ~ 1                         - diffuse 
% 3. p(sigmavc)           ~ sigmavc^-2                - diffuse
% 4. l2                   ~ iG(v,s)                   - see Gelman et al (2003), Gelman (2006)
% 5. L2c is non-stochastic and calculated from the data.

% Jarocinski (2010) follows Gelman 2006 and uses v=-1 s=0, equivalent to a
% prior of U[0,infinity] for the standard deviations. This is an improper
% distribution (as the DF are negative) but leads to a proper posterior
% below. 

% I set vbar=M+3, there isn't a permissable diffuse prior for this
% hyperparameter that conjugates. This number guarentees the prior mean
% covariance exists while imposing minimal shrinkage.
%
%---------------------------- Identification ------------------------------
% Identification is conducted via a narrative instrument such that:
% Ic=Uc*upsilonc+Vc
% upsilonc can be used to identify the structural shocks.     

% In the benchmark specification of Bahaj (2019) Vc ~ t distribution with vu degrees of
% freedom. If vu<100 this algorithm studentises the errors associated with
% the identification stage. However, studentising the errors is
% computatonally intensive. Therefore, if vu>100 the algorthim uses a
% normal prior for Vc. The matrix transforms associated with this are
% simpler; thus the algorithm runs faster.
%
%------------------------------- Inputs -----------------------------------

% Yraw is a TxMxN matrix containing endogenous VAR data,where T is the number of time series
% observations, while M is the number of VAR variables and N is the number
% of cross-sections in the panel.

% Xraw is a TxGxN matrix containing exogenous VAR data (for a VARX model),where T is the number of time series
% observations, while G is the number of exogenous VAR variables and N is the number
% of cross-sections in the panel.

% Iraw is a Tx1xN matrix containing the proxy variable for the shock of interest. 
% specify the Iraw such that NaNs at the beginning and the end of the data
% denote observations where the proxy is not available

% Define specification of the VAR model
p=varparaminputs.p;         % Number of lags on dependent variables
q=varparaminputs.q;         % Lags on exogenous variables to include, 2x1: e.g [0,2] implies contemporaneous and two lags - can't use leads.
detr=varparaminputs.detr;   % 0: no determinstic, 1: intercepts, 
                            % 2: intercept+linear trend, 3:intercept+quad trend
vu=varparaminputs.vu;       % degrees of freedom on identification t shocks
ihor=varparaminputs.ihor;   % Horizon to compute impulse responses
nburn=varparaminputs.nburn; % size of the burn-in the Gibbs-Sampler
thin=varparaminputs.thin;   % thinning factor, only save one in every thin draws
nsave=varparaminputs.nsave; % number of saved draws from the posterior
initial=varparaminputs.initial; %1= using country LS, 0=using pooled LS for initial values 

%------------------------- definition of outputs --------------------------

% These are stored in the structure "outputs" 

% Sampled posteriors
% betac_draws     =  a sample of draws of betac
% betabar_draws   =  a sample of draws of betabar 
% gammac_draws    =  a sample of draws of gammac 
% SIGMAc_draws    =  a sample of draws of SIGMAc 
% SIGMA_draws     =  a sample of draws of SIGMA 
% l1_draws        =  a sample of draws of l1      
% Uc_draws        =  a sample of draws of the residuals
% ac_draws        =  a sample of draws of the a1 identification vector
% abar_draws      =  a sample of draws of the a1bar identification vector
% bc_draws        =  a sample of draws of the b1 identification vector
% bcbar_draws     =  a sample of draws of the b1bar identification vector
% ec_draws        =  a sample of draws of the identified shocks
% l2_draws        =  a sample of draws of l2  
% relyc_draws     =  a sample of draws of the country reliability estimator
% relybar_draws   =  a sample of draws of the average reliability estimator
% Vc_draws        =  a sample of draws of identification residuals 
% Xic_draws       =  a sample of draws of t-shocks

% IRFS:
% sirs         = sorted impulse responses. (M,1,ihor,numdraws)
% Note that impulses are calculated using betabar and rotations of a pooled
% across country covariance matrix. The rows denote the shocked variable.

%----------------------------- Sampling ----------------------------------- 
madness = 0.2*nburn;        % Draws before sampler starts filtering explosive draws - switch this off by setting madness=nsave+nburn
ntot = nsave*thin + nburn;  % Total number of draws 
it_print = 100;             % Print on the screen every "it_print"-th iteration

%------------------------- Data Handling ----------------------------------
% Get initial dimensions of dependent variable
[Traw M N] = size(Yraw);  % length of raw time series, number of endogenous series, number of countries
G=size(Xraw,2);           % number of exogenous series - maybe zero
k=max(p,max(q));          % maximum lag order

% Traw was the dimension of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - k;
% Total number of observations
Obs=T*N;

%Stacking
[Y,X] = PVARX_str_H(Yraw,Xraw,detr,p,q);
I     = Iraw((k+1):end,:,:);

% Get the number of parameters per equation
[K] = size(X,2);
% number of exogenous slope parameters per equation
SX=(max(q)-min(q)+1)*G;
% number of dependent slope parameters per equation
SY=p*M;
% Number of deterministic parameters per equation
D=K-SY-SX;
% Total number of parameters
NP=M*K;   
% Total number of deterministic parameters
ND=M*D;   % this will be the size of gamma
% Total number of slope parameters
NB=NP-ND;   % this will be the size of beta

Xc=X(:,(D+1):end,:); %Just slope terms
Zc=X(:,1:D,:);       %Just deterministic terms


% ---------------- Get length of the proxy variable availability ----------
if max(isnan(I(:,:,1)))==0
t1=1;
t2=T;
Tm=t2-t1+1;
elseif isnan(I(1,:,1))==0&&isnan(I(end,:,1))~=0 
t1=1;
t2=find(diff(isnan(I(:,:,1))));
Tm=t2-t1+1;
elseif isnan(I(end,:,1))==0&&isnan(I(1,:,1))~=0 
t1=find(diff(isnan(I(:,:,1))))+1;    
t2=T;   
Tm=t2-t1+1;
else
t1=find(diff(isnan(I(:,:,1))));
t2=t1(2);
t1=t1(1)+1;
Tm=t2-t1+1;
end
% Notes: t1 is the firsy observation where the proxy is nonmissing, t2 is
% the last observation. The code can only handle missing observations at
% the start or end of the sample not gaps in proxy availability. Tm is the
% total number of time series observations where the time series observed.

%-------------------- Country by Country OLS Estimation -------------------
ALPHAc_OLS   = zeros(K,M,N);
BETAc_OLS    = zeros(K-D,M,N);
GAMMAc_OLS   = zeros(D,M,N); 
SIGMAc_OLS   = zeros(M,M,N);
Uc_OLS       = zeros(size(Y));
Ic           = zeros(Tm,1,N);
UPSILONc_OLS = zeros(M,1,N);
sigmavc_OLS  = zeros(1,1,N);
Vc_OLS       = zeros(size(Ic));
ac_OLS       = zeros(1,M,N);
bc_OLS       = zeros(M,1,N);

% First get LS estimators
for i=1:N
ALPHAc_OLS(:,:,i)   =  inv(X(:,:,i)'*X(:,:,i))*(X(:,:,i)'*Y(:,:,i)); % This is the matrix of regression coefficients, as an array for each county
BETAc_OLS(:,:,i)    =  ALPHAc_OLS((D+1):end,:,i);                    % This is the matrix of slope coefficients, as an array for each county
GAMMAc_OLS(:,:,i)   =  ALPHAc_OLS(1:D,:,i);                          % This is the matrix of deterministic coefficients, as an array for each county
Uc_OLS(:,:,i)       =  (Y(:,:,i) - X(:,:,i)*ALPHAc_OLS(:,:,i));      % Reduced Form Residuals
SSE                 =  Uc_OLS(:,:,i)'*Uc_OLS(:,:,i);                 % Sum of squared errors
SIGMAc_OLS(:,:,i)   =  SSE./(T-K+1);                                 % This is the covariance matrix as an array for each county       
Ic(:,:,i)           =  I(t1:t2,:,i);                                 % The proxy for the identification sample period   
%Ic(:,:,i)  =I(t1:t2,:,i)-X(t1:t2,:,i)*inv(X(t1:t2,:,i)'*X(t1:t2,:,i))*X(t1:t2,:,i)'*I(t1:t2,:,i);  %Project instrument on VAR DATA (comment out to switch off)
%Ic(:,:,i)=I(t1:t2,:,i)-mean(I(t1:t2,:,i),1);                                                       %Take out the instrument means (comment out to switch off)
UPSILONc_OLS(:,:,i) =  inv(Uc_OLS(t1:t2,:,i)'*Uc_OLS(t1:t2,:,i))*(Uc_OLS(t1:t2,:,i)'*Ic(:,:,i));    %identification regression
Vc_OLS(:,:,i)       =  Ic(:,:,i) - Uc_OLS(t1:t2,:,i)*UPSILONc_OLS(:,:,i);                           % identification error
sigmavc_OLS(:,:,i)  =  Vc_OLS(:,:,i)'*Vc_OLS(:,:,i)./(Tm-M);                                        % indentification error variance                
scale               =  inv(SIGMAc_OLS(:,:,i)); scale=scale(1,1)^(1/2)/UPSILONc_OLS(1,1,i);          %scaling parameter
ac_OLS(:,:,i)       =  scale*UPSILONc_OLS(:,:,i)';
bc_OLS(:,:,i)       =  SIGMAc_OLS(:,:,i)'*ac_OLS(:,:,i)';
end
betac_OLS  = (reshape(BETAc_OLS,NB,1,N));      % This is the vectorised version of BETA, as an array for each county
gammac_OLS = (reshape(GAMMAc_OLS,ND,1,N));     % This is the vectorised version of GAMMA, as an array for each county

% Pooled Variance estimator
SIGMA_OLS=zeros(M,M);
for i=1:N
SIGMA_OLS=SIGMA_OLS+(1/Obs)*(T)*SIGMAc_OLS(:,:,i);
end   


%---------------------------- Pooled OLS Estimation -----------------------

[Y_p,X_p] = PVARX_str(Yraw,Xraw,detr,[1,1,1],p,q);
GAMMA_p   = zeros(D*N,M); 
SIGMA_p   = zeros(M,M);
UPSILON_p = zeros(M,1);

% reduced form
ALPHA_p    =  inv(X_p'*X_p)*(X_p'*Y_p);
GAMMA_p    =  ALPHA_p(1:N*D,:);
BETA_p     =  ALPHA_p((N*D+1):end,:);
beta_p     =  reshape(BETA_p,NB,1);
U_p        =  (Y_p - X_p*ALPHA_p);
SIGMA_p    =  U_p'*U_p./(T*N-D*N+SX+SY+1);

% indicator to delete missing obs
cut=zeros(T,1);
cut([1:T]<t1,:)=1;
cut([1:T]>t2,:)=1;
cut=kron(ones(N,1),cut);
I_p           = I(:);
I_p(cut==1,:) = [];
Um_p          = U_p;
Um_p(cut==1,:)= [];   

% identifications step
UPSILON_p  =  inv(Um_p'*Um_p)*(Um_p'*I_p);
V_p        =  I_p - Um_p*UPSILON_p;      
sigmav_p   =  V_p'*V_p./(Tm*N-M);
scale      =  inv(SIGMA_p); scale=scale(1,1)^(1/2)/UPSILON_p(1,1);
a_p        =  scale*UPSILON_p(:,:)';
b_p        =  SIGMA_p*a_p';
% First get pooled LS estimators

%-------------------- Initialize the Gibbs sampler ------------------------
% Use LS values where possible for VAR parameters
if initial==0 %initial=0 implies pooled LS initialisation
for i=1:N
betac(:,:,i)     = beta_p;     
BETAc(:,:,i)     = BETA_p;
SIGMAc(:,:,i)    = SIGMA_p;    
UPSILONc(:,:,i)  = UPSILON_p;
sigmavc(:,:,i)   = sigmav_p;
Vc(:,:,i) = V_p((1+(i-1)*T:i*T),1);
ac(:,:,i) = a_p;
bc(:,:,i) = b_p;

g=[];
for j=1:detr
g=[g;GAMMA_p(i+(j-1)*M,:)];
end
GAMMAc(:,:,i)    = g;  
end
gammac    = (reshape(GAMMAc,ND,1,N));
SIGMA     = SIGMA_p;
Sbar      = T*SIGMA;
betabar   = beta_p;
UPSILONbar= UPSILON_p;
else      %initial=0 implies country by country LS initialisation
betac     = betac_OLS;     
BETAc     = BETAc_OLS;
gammac    = gammac_OLS;
GAMMAc    = GAMMAc_OLS;  
SIGMAc    = SIGMAc_OLS;    
SIGMA     = SIGMA_OLS;
Sbar      = T*SIGMA;
betabar   = mean(betac,3);
UPSILONc  = UPSILONc_OLS;
UPSILONbar= mean(UPSILONc_OLS,3);
sigmavc   = sigmavc_OLS;
Vc        = Vc_OLS;
ac        = ac_OLS;
bc        = bc_OLS;           
end


% Parameter covariance matix hyperparamters
Xic     =zeros(Tm,Tm,N);              % Initial t-shocks
l1      = 0.0001;                   % Initial  shrinkage parameter of reduced form slopes
L1c     = Lc_genX(Yraw,Xraw,p,q);   % Relative variance matrix generated from AR residuals
invL1c   = zeros(size(L1c)); 
l2      = 0.0001;                   % Initial  shrinkage parameter of reduced form slopes
L2c     = zeros(M,M,N); 
invL2c   = zeros(M,M,N); 
for i=1:N
Xic(:,:,i)=eye(Tm);
invL1c(:,:,i)=inv(L1c(:,:,i));                            % Inverted Relative variance parameter 
L2c(:,:,i)=inv(diag(var(Uc_OLS(:,:,i))/var(Ic(:,:,i))));  % Relative variance parameter 
invL2c(:,:,i)=(diag(var(Uc_OLS(:,:,i))/var(Ic(:,:,i))));  % Inverted Relative variance parameter 
end    
output.L1c=L1c;
output.L2c=L2c;

% checking stability of LS estimators - we only care about this if some
% models are explosive (can be ignored, for debugging only).
LS_stab=[];
for i=1:N
LS_stab=[LS_stab;max(abs(eig(companion(BETAc(:,:,i),p,M))))];
end


% Fixed prior hyperparameters:
s=0.00000;   % Scale parameter for the IG distribution of l 
v=-1;        % Shape parameter for the IG distribution of l
vbar=M+3;    % DF parameter for the IW distribution of SIGMAc   
%-------------------------- Storage Matrices ------------------------------

% Storage space for posterior draws
output.betac_draws   =  zeros(NB,1,N,nsave);   % save draws of betac
output.betabar_draws =  zeros(NB,1,nsave);     % save draws of betabar
output.gammac_draws  =  zeros(ND,1,N,nsave);   % save draws of gammac
output.SIGMAc_draws  =  zeros(M,M,N,nsave);    % save draws of SIGMAc
output.SIGMA_draws   =  zeros(M,M,nsave);      % save draws of SIGMAc  
output.l1_draws      =  zeros(nsave,1);        % save draws of l1 
output.Uc_draws      =  zeros(size(Uc_OLS,1),size(Uc_OLS,2),N,nsave); % save draws of residuals 
output.ac_draws      =  zeros(1,M,N,nsave);    % save draws of ac 
output.abar_draws    =  zeros(1,M,nsave);      % save draws of abar
output.UPSILONc_draws      =  zeros(1,M,N,nsave);    % save draws of UPSILONc 
output.UPSILONbar_draws    =  zeros(1,M,nsave);      % save draws of UPSILONbar 
output.bc_draws      =  zeros(M,1,N,nsave);    % save draws of bc 
output.bbar_draws    =  zeros(M,1,nsave);      % save draws of bbar 
output.ec_draws      =  zeros(T,1,N,nsave);    % save draws of identified shock
output.sigmavc_draws =  zeros(N,nsave);        % save draws of identification step se
output.l2_draws      =  zeros(nsave,1);        % save draws of l2 
output.relyc_draws   =  zeros(nsave,N);        % save draws of reliability 
output.Vc_draws      =  zeros(size(Vc_OLS,1),size(Vc_OLS,2),N,nsave); % save draws of residuals 
output.Xic_draws     =  zeros(Tm,Tm,N,nsave);  % save draws of t-shocks


%----------------------- IMPULSE RESPONSES SET-UP:  ----------------------- 
% Create matrices to store impulses

imps      = zeros(M,1,ihor,nsave); % complete set of impulse responses (mean)
impsc     = zeros(M,1,ihor,N,nsave); % complete set of impulse responses (country)
output.decomp_draws   = zeros(M,ihor,nsave);   % decomposition storage matrix (mean)
output.decompc_draws  = zeros(M,ihor,N,nsave); % decomposition storage matrix (country)


%-------------------------- Preliminaries End -----------------------------

%% ========================= Start Sampling ===============================

tic;
disp('Number of iterations');
for irep = 1:ntot  %Start the Gibbs "loop"
    accept=0;
    complete=0;
             
    if mod(irep,it_print) == 0 % print iterations
        disp(irep);
        toc;
    end
    

%------------------------------- Estimation -------------------------------
    
    
    while not(complete==1) %Everything needs to work for the draw to be accepted (i.e. all madness filters) 
    dr=0;           
    Uc=zeros(T,M,N);   % residuals of the country models (storage matrix)

    % Building the PHI matrix of joint covariances (Preliminary step).
%   if vu<100  % If the errors are studentised we need a large matrix
    PHI11=zeros(M*T,M*T,N);
    PHI12=zeros(M*T,Tm,N);
    PHI21=zeros(Tm,M*T,N);
    PHI22=zeros(Tm,Tm,N);
    SIGMAconinv=zeros(size(PHI11));
    ycon=zeros(M*T,1,N);
    
    for i=1:N
    PHI11(:,:,i)=kron(SIGMAc(:,:,i),eye(T));
    PHI12(:,:,i)=kron(SIGMAc(:,:,i)*UPSILONc(:,:,i),[zeros(t1-1,Tm);eye(Tm);zeros(T-t2,Tm)]);
    PHI21(:,:,i)=PHI12(:,:,i)';
    PHI22(:,:,i)=kron(UPSILONc(:,:,i)'*SIGMAc(:,:,i)*UPSILONc(:,:,i),eye(Tm))+sigmavc(:,:,i)*Xic(:,:,i);
    
    SIGMAconinv(:,:,i)=inv(PHI11(:,:,i)-PHI12(:,:,i)/PHI22(:,:,i)*PHI21(:,:,i)); % inverse variance of the data conditional on the proxy 
    ycon(:,:,i)=PHI12(:,:,i)/PHI22(:,:,i)*Ic(:,:,i);  % Expectation of the data conditional on the proxy
    end
    
   
    % PART A: The first level parameters in the reduced form VAR    
    
    % A1) Draw betac, betac|parameters,Data ~ Normal((Dc)^(-1)*dc,inv(Dc)) 
    stable=0;
    stable_N=zeros(N,1);
    while not(stable==1)
    for i=1:N
    stable_N(i,:)=0;
    stacker=(Y(:,:,i)-Zc(:,:,i)*GAMMAc(:,:,i));
    
    % Generate conditional density parameters
    if vu<100
    Dc=kron(eye(M),Xc(:,:,i))'*SIGMAconinv(:,:,i)*kron(eye(M),Xc(:,:,i))+(1/l1)*invL1c(:,:,i);    
    dc=kron(eye(M),Xc(:,:,i))'*SIGMAconinv(:,:,i)*(stacker(:)-ycon(:,:,i))+(1/l1)*invL1c(:,:,i)*betabar;
    else
    Dc=kron(SIGMAconinv(:,:,i),Xc(:,:,i)'*Xc(:,:,i))+(1/l1)*invL1c(:,:,i);
    dc=kron(SIGMAconinv(:,:,i),Xc(:,:,i))'*(stacker(:)-ycon(:,:,i))+(1/l1)*invL1c(:,:,i)*betabar;  
    end
      
    % Draw the slope parameters
    while not(stable_N(i,:)==1) %keep drawing until stability test passes.
    betac(:,:,i)=Dc\dc+chol(inv(Dc))'*randn(NB,1); % Draw of betac
    BETAc(:,:,i)=reshape(betac(:,:,i),K-D,M);           % rehaping into matrix
    % stability test
    if irep<madness
    stable_N(i,:)=1;
    else
    BB=companion(BETAc(:,:,i),p,M);                     
    stable_N(i,:)=not(max(abs(eig(BB)))>1.00);
    end
    %if stable==0
    %disp(['country ', num2str(i),' is unstable']);    
    %end
    
    end
    stable_N(i,:)=0;
    end
    stable=min(stable_N); 
      
    % A2) Draw betabar, betabar|parameters,Data ~ Normal((Gc)^(-1)*gc,inv(Gc)) 
    eigs=[];
    Gc=zeros(NB,NB);
    gc=zeros(NB,1);
    for i=1:N
    Gc=Gc+(1/l1)*invL1c(:,:,i);    
    gc=gc+(1/l1)*invL1c(:,:,i)*betac(:,:,i); 
    end
    Gbar=inv(Gc)*gc;
    GG=companion(reshape(Gbar,K-D,M),p,M);  
    
    if (max(abs(eig(GG)))<1) % mean is stable
    %keep drawing until stability test passes.    
    while not(stable==1)
    betabar=Gbar+chol(inv(Gc))'*randn(NB,1); % Draw of betabar
    BETAbar=reshape(betabar,K-D,M);          % rehaping into matrix
    % stability test
    if irep <madness
    stable=1;
    else
    BB=companion(BETAbar,p,M);                     
    stable=not(max(abs(eig(BB)))>1);
    if stable==0
    eigs=[eigs;max(abs(eig(BB)))];
    %disp('drawn common mean is not stable');
    end    
    end
    end
    
    else
    %disp('distribution mean is unstable')
    stable=0;
    end
    
    end
    
    
    % A3) Draw gammac, gammac|parameters,Data ~ Normal((Fc)^(-1)*fc,inv(Fc)) 
    
    for i=1:N
    stacker=Y(:,:,i)-Xc(:,:,i)*BETAc(:,:,i);    
    if vu<100
    Fc=kron(eye(M),Zc(:,:,i))'*SIGMAconinv(:,:,i)*kron(eye(M),Zc(:,:,i));
    fc=kron(eye(M),Zc(:,:,i))'*SIGMAconinv(:,:,i)*(stacker(:)-ycon(:,:,i));
    else
    Fc=kron(SIGMAconinv(:,:,i),Zc(:,:,i)'*Zc(:,:,i));
    fc=kron(SIGMAconinv(:,:,i),Zc(:,:,i))'*(stacker(:)-ycon(:,:,i));      
    end
    gammac(:,:,i)=inv(Fc)*fc+chol(inv(Fc))'*randn(ND,1); % Draw of gammac
    GAMMAc(:,:,i)=reshape(gammac(:,:,i),D,M);            % rehaping into matrix   
    end 

   
    
    % A4) Draw SIGMAc |parameters,Data ~ iW((Uc)'*(Uc)+S,T+vbar) 
    for i=1:N
    Uc(:,:,i)=Y(:,:,i)-Xc(:,:,i)*BETAc(:,:,i)-Zc(:,:,i)*GAMMAc(:,:,i);  
    SIGMAc(:,:,i) = iwishrnd((Uc(:,:,i)'*Uc(:,:,i))+Sbar,T+vbar);% Draw SIGMA    
    end    
    
    
     % A5) Draw S |parameters,Data ~ iW(Sc,N*vbar)
    Sc=zeros(M,M);
    for i=1:N
    Sc=Sc+inv(SIGMAc(:,:,i));
    end    
    Sbar = wishrnd(inv(Sc),N*vbar);% Draw SIGMA  
    
    SIGMA=Sbar/vbar;
    
    % PART B: Draw the hyperparameters for the reduced form VAR
    % B1) Draw l1, l1|parameters,Data ~ iG(s+Kc,v+M*N*P) 
    Kc=0;
    for i=1:N
    Kc=Kc+(betac(:,:,i)-betabar)'/(L1c(:,:,i))*(betac(:,:,i)-betabar);    
    end    
    l1=1/gamrnd((0.5*(v+N*M*M*p)),(0.5*(s+Kc))^(-1));
    
    % PART C: Draw the parameters for the for the indetification model:    
    
    % C1) Draw UPSILONc, UPSILONc|parameters,Data ~ Normal((Jc)^(-1)*jc,inv(Jc)) 
    for i=1:N
    Uc(:,:,i)=Y(:,:,i)-Xc(:,:,i)*BETAc(:,:,i)-Zc(:,:,i)*GAMMAc(:,:,i);  
    end
    
    if vu<100
    for i=1:N
    Ucc=Uc(t1:t2,:,i);   % Keep means
    Jc=1/sigmavc(:,:,i)*Ucc'/Xic(:,:,i)*Ucc+(1/l2)*invL2c(:,:,i);    
    stacker=Ic(:,:,i);  % Keep means
    jc=(1/sigmavc(:,:,i))*Ucc'/Xic(:,:,i)*stacker(:)+(1/l2)*invL2c(:,:,i)*UPSILONbar;
    UPSILONc(:,:,i)=Jc\jc+chol(inv(Jc))'*randn(M,1); % Draw of UPSILONc
    end
    else
    for i=1:N
    Ucc=Uc(t1:t2,:,i);
    Jc=(1/sigmavc(:,:,i))*Ucc'*Ucc+(1/l2)*invL2c(:,:,i);    
    %stacker=Ic(:,:,i)-mean(Ic(:,:,i),1);
    stacker=Ic(:,:,i);
    jc=(1/sigmavc(:,:,i))*Ucc'*stacker(:)+(1/l2)*invL2c(:,:,i)*UPSILONbar;
    UPSILONc(:,:,i)=Jc\jc+chol(inv(Jc))'*randn(M,1); % Draw of UPSILONc
    end
    end

    % C2) Draw UPSILONbar, UPSILONbar|parameters,Data ~ Normal((Hc)^(-1)*hc,inv(Hc))      
    Hc=zeros(M,M);
    hc=zeros(M,1);
    for i=1:N
    Hc=Hc+(1/l2)*invL2c(:,:,i);    
    hc=hc+(1/l2)*invL2c(:,:,i)*UPSILONc(:,:,i); 
    end
    UPSILONbar=Hc\hc+chol(inv(Hc))'*randn(M,1); % Draw of UPSILONc
    
    % C3) Draw sigmavc |parameters,Data ~ iG((Uc)'*(Uc)+S,T+vbar) 
    if vu<100
    for i=1:N    
    Ucc=Uc(t1:t2,:,i);    
    Vc(:,:,i)=Ic(:,:,i)-Ucc*UPSILONc(:,:,i);  
    sigmavc(:,:,i) = iwishrnd((Vc(:,:,i)'/Xic(:,:,i)*Vc(:,:,i)),Tm);% Draw SIGMA        
    end
    else
    for i=1:N
    Ucc=Uc(t1:t2,:,i);
    Vc(:,:,i)=Ic(:,:,i)-Ucc*UPSILONc(:,:,i);  
    sigmavc(:,:,i) = iwishrnd((Vc(:,:,i)'*Vc(:,:,i)),Tm);% Draw SIGMA    
    end
    end

    
    % PART D: Draw the hyperparameters for the for the identification model:
    % D1) Draw l2, l2|parameters,Data ~ iG(s+Mc,v+M*N) 
    Mc=0;
    for i=1:N
    Mc=Mc+(UPSILONc(:,:,i)-UPSILONbar)'/(L2c(:,:,i))*(UPSILONc(:,:,i)-UPSILONbar);    
    end    
    l2=1/gamrnd((0.5*(v+N*M)),(0.5*(s+Mc))^(-1));

    % D2) Draw Xic, (see notes for distrn)
    if vu<100
    for i=1:N
    draw=random('chi2',vu+1,[Tm,1]);
    for ii=1:Tm
    Xic(ii,ii,i)=(sigmavc(:,:,i)^(-1)*Vc(ii,:,i)^2+vu)/draw(ii);
    end
    end 
    end
    
    
    % E:  Prep for impulses
    for i=1:N  
    scale           =  (UPSILONc(:,:,i)'*SIGMAc(:,:,i)*UPSILONc(:,:,i))^(1/2); 
    ac(:,:,i)       =  (1/scale)*UPSILONc(:,:,i)';
    bc(:,:,i)       =  SIGMAc(:,:,i)'*ac(:,:,i)';    
    end    
    
    scale           =  (UPSILONbar'*SIGMA*UPSILONbar)^(1/2); 
    abar            =  (1/scale)*UPSILONbar(:,:)';
    bbar            =  SIGMA(:,:)*abar(:,:)';   
    complete=1; 
    end    
    % =============Estimation ends here

    
    
%--------------------------- IMPULSE RESPONSES ----------------------------
    if irep > nburn
  
    if mod(irep-nburn,thin)==0
    draw=(irep-nburn)/thin;    
        J=zeros(M,M*p);
        J(1:M,1:M)=eye(M);
        
        % Average country model
        BB=companion(BETAbar,p,M);    
        MSE=zeros(M,M,ihor);
        decomp=zeros(M,ihor);
        impu=zeros(size(BB,1),1,ihor);
        for j=1:ihor
        % impulse responses
        if j==1
            impu(1:M,:,j)=bbar;
        elseif j>1
            impu(:,:,j)=BB*impu(:,:,j-1);
        end
        imps(:,:,:,draw)=impu(1:M,1,:);
       
        % MSE error
        theta=J*BB^(j-1)*J';
        if j==1    
        MSE(:,:,j)=theta*SIGMA*theta';
        else 
        MSE(:,:,j)=theta*SIGMA*theta'+MSE(:,:,j-1);
        end 
        
        % Var decomp
        for jj=1:M
        decomp(jj,j)=(sum(impu(jj,1,1:j).^2))/MSE(jj,jj,j);
        end
        end
       
        
        % Individual country model
        decompc=zeros(M,ihor,N);
        for i=1:N
        BB=companion(BETAc(:,:,i),p,M);    
        impu=zeros(size(BB,1),1,ihor);   
        MSE=zeros(M,M,ihor);
        
        
        for j=1:ihor
        % impulses
        if j==1
            impu(1:M,:,j)=bc(:,:,i);
        elseif j>1
            impu(:,:,j)=BB*impu(:,:,j-1);
        end
        impsc(:,:,:,i,draw)=impu(1:M,1,:);
        % MSE error
        theta=J*BB^(j-1)*J';
        if j==1    
        MSE(:,:,j)=theta*SIGMAc(:,:,i)*theta';
        else 
        MSE(:,:,j)=theta*SIGMAc(:,:,i)*theta'+MSE(:,:,j-1);
        end
        
        
        % Var decomp
        for jj=1:M
        decompc(jj,j,i)=(sum(impu(jj,1,1:j).^2))/MSE(jj,jj,j);
        end
        
        end
        end

%--------------------------- Extract Shocks -------------------------------        
ec    =  zeros(T,1,N);  
for i=1:N
ec(:,:,i)=Uc(:,:,i)*ac(:,:,i)';    
end


%--------------------- Save draws of the parameters -----------------------
        output.gammac_draws(:,:,:,draw)         = gammac;
        output.betac_draws(:,:,:,draw)          = betac;
        output.betabar_draws(:,:,draw)          = betabar;
        output.SIGMAc_draws(:,:,:,draw)         = SIGMAc;
        output.SIGMA_draws(:,:,draw)            = SIGMA;
        output.l1_draws(draw)                   = l1;
        output.Uc_draws(:,:,:,draw)             = Uc;
        output.ac_draws(:,:,:,draw)             = ac;  
        output.abar_draws(:,:,draw)             = abar;
        output.UPSILONc_draws(:,:,:,draw)       = UPSILONc;  
        output.UPSILONbar_draws(:,:,draw)       = UPSILONbar;   
        output.bc_draws(:,:,:,draw)             = bc;    
        output.bbar_draws(:,:,draw)             = bbar;      
        output.ec_draws(:,:,:,draw)             = ec;    
        output.l2_draws(draw)                   = l2;        
        output.Vc_draws(:,:,:,draw)             = Vc;
        output.Xic_draws(:,:,:,draw)            = Xic;
        output.sigmavc_draws(:,draw)            = squeeze(sigmavc);
        output.decomp_draws(:,:,draw)           = decomp;
        output.decompc_draws(:,:,:,draw)        = decompc;
        for i=1:N
        output.relyc_draws(draw,i)              = UPSILONc(:,:,i)'*SIGMAc(:,:,i)*UPSILONc(:,:,i)/(UPSILONc(:,:,i)'*SIGMAc(:,:,i)*UPSILONc(:,:,i)+sigmavc(:,:,i));         
        end
        output.relybar_draws(draw)              = UPSILONbar'*SIGMA*UPSILONbar/(UPSILONbar'*SIGMA*UPSILONbar+mean(sigmavc,3));
        
    end

    end % end saving results
       
end %end the main Gibbs for loop
%====================== End Sampling Posteriors ===========================






%======================== Impulse Responses ===============================
output.sir=sort(imps,4);
output.sirc=zeros(size(impsc));
for i=1:N
output.sirc(:,:,:,i,:)=sort(impsc(:,:,:,i,:),5);
end

end

