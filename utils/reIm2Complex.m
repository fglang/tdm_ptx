function out = reIm2Complex(x)
    x = vec(x);
    L = length(x);
    out = x(1:L/2) + 1j .* x(L/2+1:end);
end