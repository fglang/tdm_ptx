function out = complex2ReIm(x)
    x = vec(x);
    out = cat(1,real(x),imag(x));
end