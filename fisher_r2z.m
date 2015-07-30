function z=fisher_r2z(r,n)

prevst=warning('off','MATLAB:divideByZero');
%z=real(sqrt(n-3)/2 .* log((1+r)./(1-r)));
z=real(1/2 .* log((1+r)./(1-r)));
% z=sqrt(n-3)/2 * log((1+y)./(1-y));
warning(prevst);