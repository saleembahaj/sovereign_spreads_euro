function [y_counter] = counter_fac(yraw,betac,bc,ec,p,M,detr,start)


Tn=size(yraw,1);
N=size(yraw,3);
numdraws=size(betac,4);
y_counter=zeros(Tn,M,N,numdraws);

for cc=1:N
parfor i=1:numdraws
res2=[];
y=yraw(:,:,cc);    
beta=squeeze(betac(:,:,cc,i));
BETA=companion(beta,p,M);

% Lets get the regression figures
[~,X] = PVAR_str_H(y,detr,p);
% Get the number of parameters per equation
[K] = size(X,2);
% Number of deterministic parameters per equation
D=K-p*M;
% Total number of parameters
NP=M*K;   
% Total number of deterministic parameters
ND=M*D;   % this will be the size of gamma


% Generate change 
shocks=ec(:,:,cc,i);
b=bc(:,:,cc,i);
res2=shocks*b';
res2=[zeros(p,M);res2];             % adding zero time t-p observation
res2(1:start,:)=zeros(start,M);     % setting all shocks prior to the counterfactual start date to zero


% Companion forms
rcomp=zeros(Tn-p,M*p);              % new residual companion forms
ycomp=zeros(Tn-p,M*p);              % yraw companion form
for t=1:Tn-p
for j=1:p
rcomp(t,(j-1)*M+1:j*M)=res2(t+p-j+1,:);   
ycomp(t,(j-1)*M+1:j*M)=y(t+p-j+1,:);   
end
end


% Generating counterfactual series 
shock_imp=zeros(Tn-p,M*p); %impact of the shocks to be removed
for t=1:(Tn-p)    
for j=1:t    
shock_imp(t,:)=shock_imp(t,:)+((BETA^j)*rcomp(t-j+1,:)')';
end    
end
shock_imp=[zeros(p,M*p);shock_imp];
y_counter(:,:,cc,i)=y-shock_imp(:,1:M);  %And the counterfactual!!

end


end
end

