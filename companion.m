% Construct the companion representation of a state space model
% with state vector theta p lags in the VAR representation and 
% n variables.

function C = companion(theta,p,n)
if size(theta,2)<n %f theta is vectorised convert back into matrix form
theta= reshape(theta,size(theta,1)/n,n);
end
theta = theta((size(theta,1)-n*p+1):end,:);


C = [theta' ; eye(n*(p-1)) zeros(n*(p-1),n)];