function [ w, r, x, P ] = capping( w, r, x, P, threshold )

if length(w) > threshold
    [~,idxsort]= sort(w,'descend');
    idxkeep=idxsort(1:threshold);
    w= w(idxkeep);
    r = r(idxkeep);
    x = x(idxkeep);
    P = P(idxkeep);
end
w = normalize(w);

end

