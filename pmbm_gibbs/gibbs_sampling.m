function [ assignments ] = gibbs_sampling( wupd,wmiss,wnew,M )

m = length(wnew);
n = length(wmiss);
w = wupd./repmat(wmiss,1,m);
num_iteration = 100;
assignments = zeros(num_iteration,m);
%Initialize with measurements are all assigned to new tracks
%Loop through measurements
unused_tracks = true(n,1); %indices of unused tracks
for i = 1:num_iteration
    for j = 1:m
        if i>=2 && assignments(i-1,j)>0
            unused_tracks(assignments(i-1,j)) = true; %release track
        end
        weight = [w(unused_tracks,j);wnew(j)];
        cdf = cumsum(weight/sum(weight));
        I = find(cdf>rand(),1); %sampling
        if I~=length(cdf)
            idx = find(unused_tracks==true,I);
            idx = idx(end);
            unused_tracks(idx) = false;
            assignments(i,j) = idx;
        else
            assignments(i,j) = 0;% corresponds to assignment to new track
        end
        
    end
end
if size(assignments,1) > M
    assignments = assignments(end-M+1:end,:);
end
assignments = unique(assignments,'rows');


end

