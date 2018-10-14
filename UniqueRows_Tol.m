function [T1,ia] = UniqueRows_Tol(Xin, tol)

% check for/supply defaults
x = Xin; [n,p] = size(x);

if ((nargin<2) || isempty(tol))
      tol = 1.e-12;
end

tol = tol*(1+10*eps); %%% Why?

% consolidate elements of x.
% first shift, scale, and then round up. 
xhat = x - repmat(min(x,[],1),n,1)+tol*eps;
xhat = ceil(xhat/tol);


[xhat,tags] = sortrows(xhat);

% count the replicates
iu = [true;any(diff(xhat),2)];
ia = sort(tags(iu == 1));
T1 = Xin(ia,:);