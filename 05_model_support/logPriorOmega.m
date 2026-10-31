function funval = logPriorOmega(model_specs)

Omega = model_specs.omega{1}(1);
aa = 2;
bb = 2;

if ~isscalar(Omega) || ~isreal(Omega) || ~isfinite(Omega) || ...
        Omega <= 0 || Omega >= 1
    funval = -Inf;
    return;
end

funval = (aa-1)*log(Omega) + ...
    (bb-1)*log1p(-Omega) - betaln(aa,bb);
end
