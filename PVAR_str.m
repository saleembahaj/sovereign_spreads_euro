%--------------------------------------------------------------------------
% PVAR_str
% Creates matrices for PVAR(k), estimation via LS
% 
%-------------------------------------------------------------------------- 

% Inputs:
% y: TxNxM Panel time series observations of data 
% c=0 implies no constant
% c=1 implies ficed effects and no trend
% c=2 implies fixed effects and linear trend
% c=3 implies fixed effects and quadtratic linear trends
% p is panel assumptions: e.g p=[1 0 0] implies unit specific constants and
% common slopes coefficients.
% k: lag length  

% Outputs:
% yy: set of dependent variables
% xx: set of explanatory variables such that (x'x)^(-1)*x'yy gives coeffs


function [yy,xx] = PVAR_str(y,c,p,k)
%s=size(y);
T=size(y,1); N=size(y,2); M=size(y,3);


%Lagging variables 
yy=zeros(T-k,N,M);
xx=zeros(T-k,N*k,M);
for ii=1:M
for i=1:N,
yy(:,i,ii)=y(k+1:T,i,ii);
end
for j=1:k,   
xx(:,((j-1)*(N)+1):N*j,ii)=y(k+1-j:T-j,1:N,ii);
end
end

%stacking the cross section
yy2=[];
xx2=[];
for i=1:M
yy2=[yy2;yy(:,:,i)];
xx2=[xx2;xx(:,:,i)];
end
yy=yy2;
xx=xx2;


% Creating deterministic components
if c>0
if p(:,1)==1;    
cons=kron(eye(M),ones(T-k,1));
else
cons=kron(ones(M,1),ones(T-k,1));
end
if c>1
if p(:,2)==1;    
lin=kron(eye(M),(1:T-k)');
else
lin=kron(ones(M,1),(1:T-k)');
end
if c>2
if p(:,2)==1;    
quad=kron(eye(M),((1:T-k).^2)');
else
quad=kron(ones(M,1),((1:T-k).^2)');
end
end
end
end


% Adding deterministic components
if c==1
xx=[cons xx];
elseif c==2
xx=[cons lin xx];
elseif c==3    
xx=[cons lin quad xx];       
end


    
end



% This is a modified version of the code written by:
% L.Gambetti and F. Canova, version 2.1, 6-10-2005
% Modified by Saleem Bahaj August 2012