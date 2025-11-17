function res = applyMatrix2nd(mat, nd, ax)
% apply 2d matrix mat along ax's axis of nd array
% faster and less memory consuming than einsum!

    if nargin < 3
        ax=1;
    end
    sz = size(nd);
    dims = 1:ndims(nd);

    neworder = [ax, setdiff(dims,ax)]; % permutation indices to bring axis that is used for matrix product to front

    nd = permute(nd,neworder); % ax, rest
    nd = reshape(nd,sz(ax),[]);
    res = mat*nd;
    res = reshape(res, [size(res,1), sz(neworder(2:end))]);
    [~,isort] = sort(neworder); % undo permutation indices
    res = permute(res, isort);
end