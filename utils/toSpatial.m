function out = toSpatial(inp, mask)
% convert masked list of things inp [nspat, other dims] into full spatial
% thing according to mask [nx,ny,nz] by filling spatial axis into mask
% result will be [nx,ny,nz, other dims]
% mask in NaN convention
    sz = size(mask);

    if islogical(mask)
        mask = double(mask);
        mask(mask==0) = NaN;
    end

    szin = size(inp);
    szrest = szin(2:end);
    out = NaN([prod(sz), prod(szrest)]);
    out(~isnan(mask),:) = inp;
    out = reshape(out,[sz, szrest]);
end