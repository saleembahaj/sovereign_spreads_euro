%Generates a prior covariance matrix for VAR slope parameters along the
%lines describe in Jarocinski (2010):

% For reference:
% betac ~ Normal(beta,l*Lc)
% beta and l are determined outside this function.

% Lc is such that Cov(betac(i),betac(j))=0 if i=/=j and
% Var(betac)=l*sn^2/sk^2. Where n is the variable corresponing to the
% equation being estimated and k is the variable which the coefficient is
% on. For example, the variance of the coefficient on the variables own lag
% is equal to l always. These variances are the same regardless of the lag
% under consideration. sn^2 corresponds to the variance of the residuals of
% an AR(p) evaluated on variable n.



function [Lc] = Lc_genX(Yraw,Xraw,p,q)
[Traw M N] = size(Yraw);
G=size(Xraw,2);

k=max(q);
Lc=zeros(M*(G*(k-min(q)+1)+p*M),M*(G*(k-min(q)+1)+p*M),N); % storage matrix for Lc

for i=1:N

% Creating variances of residual processes.
varar=zeros(M,1);
for j=1:M
[Y,X] = PVAR_str_H(Yraw(:,j,i),0,p);
B=inv(X'*X)*(X'*Y);
varar(j)=(Y - X*B)'*(Y - X*B)./(Traw-2*p);   %estimated variance 
end    

vararX=[];
if not(isempty(Xraw))
k=k+1-min(q);
vararX=zeros(G,1);
for j=1:G
[Y,X] = PVAR_str_H(Xraw(:,j,i),0,k);
B=inv(X'*X)*(X'*Y);
vararX(j)=(Y - X*B)'*(Y - X*B)./(Traw-2*k);   
end
end

k=k-min(q)+1;

% Generating diagonal elements
varlist=[kron(ones(k,1),vararX);kron(ones(p,1),varar)];

diagL=kron(ones(M,1)./varar,varlist); %Creating the diagnol

Lc(:,:,i)=inv(diag(diagL));


end

end

