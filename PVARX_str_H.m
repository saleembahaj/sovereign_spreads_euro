%--------------------------------------------------------------------------
% Pvar_str_H
% Creates matrices for VARX(p,q), handles Panel data but retains the cross
% section dimension rather than stacking
% 
%-------------------------------------------------------------------------- 

% Inputs:
% y: TxNxM Panel time series observations of data
% x= TxGxM Panel time series observations of exogenous data.
% c=0 implies no constant
% c=1 implies constant and no trend
% c=2 implies constant and linear trend
% c=3 implies constant and quadtratic linear trends
% p: lag length on dependent variables
% q: 2x1 vector if q=[0,n] included first n lags of x plus contemporaneous
% term, if q=[2,n] include the 2nd through nth lags of x.

% Outputs:
% yy: set of dependent variables
% xx: set of explanatory variables such that (x'x)^(-1)*x'yy gives coeffs


function [yy,xx] = PVARX_str_H(y,x,c,p,q)
if isempty (x)
[yy,xx] = PVAR_str_H(y,c,p);    
else
%s=size(y);
T=size(y,1); N=size(y,2); M=size(y,3); G=size(x,2);

if not(size(x,3)==size(y,3));
error('exogenous variables do not have the same cross section dimensions as dependent variables') 
end    

k=max(p,max(q));

% Creating deterministic components
det=[];
if c==1   
det=ones(T-k,1);
elseif c==2   
det=[ones(T-k,1),(1:T-k)'];
elseif c==3
det=[ones(T-k,1),(1:T-k)',((1:T-k).^2)'];    
end



%Lagging variables 
yy=zeros(T-k,N,M);
ex=zeros(T-k,G*(max(q)-min(q)+1),M);
xx=zeros(T-k,N*k,M);
xx2=[];
for ii=1:M
for i=1:N,
yy(:,i,ii)=y(k+1:T,i,ii);
end
for j=1:p,   
xx(:,((j-1)*(N)+1):N*j,ii)=y(k+1-j:T-j,1:N,ii);
end
for j=min(q):1:max(q)
jj=j;
    if min(q)==0
    jj=j+1;
    end
    ex(:,((jj-1)*(G)+1):G*jj,ii)=x(k+1-j:T-j,1:G,ii);
end    
xx2=cat(3,xx2,[det,ex(:,:,ii),xx(:,:,ii)]);
end
xx=xx2;

end    
end



% This is a modified version of the code written by:
% L.Gambetti and F. Canova, version 2.1, 6-10-2005
% Modified by Saleem Bahaj October 2012