function [Y X0] = filt_tc(y,TR,cut_off,linear)

if nargin<4
    linear = 0;
end
    
    
[T k] = size(y);

n       = fix(2*(T*TR)/cut_off + 1);
X0      = spm_dctmtx(T,n);
X0      = X0(:,2:end);

if linear, % linear trend
    X0 = [linspace(max(X0(:)),min(X0(:)),size(X0,1))', X0];
end

for i = 1:k
beta    = X0\y(:,i);

Y(:,i)  = y(:,i) - X0*beta;
end