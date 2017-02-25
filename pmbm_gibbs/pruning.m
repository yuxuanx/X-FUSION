function [ w, r, x, p ] = pruning( w, r, x, p, threshold )

if ~isempty(w)
    idxPrune = w > threshold;
    w = w(idxPrune);
    r = r(idxPrune);
    x = x(idxPrune);
    p = p(idxPrune);
    w = normalize(w);
end

end

