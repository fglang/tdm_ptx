function out = nrmse2(a,b)
% normalize by L2 norm of a, such that nrmse=L2(a-b)/L2(a)
% -> not commutative! nrmse(a,b)~=nrmse(b,a) !
L2 = @(x) sqrt(mean(abs(x).^2,"all","omitnan"));
out = L2(a-b) ./ L2(a);