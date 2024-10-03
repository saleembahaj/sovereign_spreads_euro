%--------------------------------------------------------------------------
% Pvar_str_H
% Creates matrices for VAR(k), handles Panel data but retains the cross
% section dimension rather than stacking
% 
%-------------------------------------------------------------------------- 

% Inputs:
% y: TxNxM Panel time series observations of data 
% c=0 implies no constant
% c=1 implies constant and no trend
% c=2 implies constant and linear trend
% c=3 implies constant and quadtratic linear trends
% k: lag length  

% Outputs:
% yy: set of dependent variables
% xx: set of explanatory variables such that (x'x)^(-1)*x'yy gives coeffs
% x: useless as far as I can see 

function [yy,xx] = PVAR_str_H(y,c,k)
%s=size(y);
T=size(y,1); N=size(y,2); M=size(y,3);


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
xx=zeros(T-k,N*k,M);
xx2=[];
for ii=1:M
for i=1:N,
yy(:,i,ii)=y(k+1:T,i,ii);
end
for j=1:k,   
xx(:,((j-1)*(N)+1):N*j,ii)=y(k+1-j:T-j,1:N,ii);
end
xx2=cat(3,xx2,[det,xx(:,:,ii)]);
end
xx=xx2;

    
end



% This is a modified version of the code written by:
% L.Gambetti and F. Canova, version 2.1, 6-10-2005
% Modified by Saleem Bahaj October 2012