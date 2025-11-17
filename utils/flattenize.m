function B1c = flattenize(x,ax)
if nargin < 2; ax=4; end
nc = size(x,ax);
B1c = reshape(x(~isnan(x)),[],nc);