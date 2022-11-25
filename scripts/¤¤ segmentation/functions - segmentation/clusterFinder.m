function [JmaxVals,maxVals,dists,C] = clusterFinder(z,u,r)
% Takes a vector of values z and a threshold u and clusters the data-points
% that lie above u. r is the number of consecutive values falling below u
% required to initiate a new cluster.
%
% Output:
% maxInd  = indeces of the maximum of each respective cluster
% maxVals = maximum values of each cluster
% dists   = distances between the peaks

uInd = find(z>u);
% vector that keeps cluster label for each threshold exceedence
clusters = zeros(1,length(uInd));
% start at cluster 1
clusters(1) = 1;
% loop variable that keeps track of current cluster value
cluster = 1;

for k=1:length(uInd)-1
    diff = uInd(k+1)-uInd(k);
    if diff < r
        clusters(k+1) = cluster;
    else
        cluster = cluster + 1;
        clusters(k+1) = cluster;
    end
        
end


clusterMaxInd = zeros(1,cluster);
clusterMaxVals = zeros(1,cluster);
for k=1:cluster
    [m,ind] = max(z(uInd(clusters==k)));
    ck_ind = uInd(clusters==k);
    clusterMaxInd(k) = ck_ind(ind);
    clusterMaxVals(k) = m;
end

minMax = ones(cluster,2);
for i=1:cluster
    xtemp = uInd(clusters==i);
    minMax(i,:) = [min(xtemp), max(xtemp)];
end

C.ind = uInd;
C.label = clusters;
C.minMax = minMax;

JmaxVals = clusterMaxInd;
maxVals = clusterMaxVals;
dists = JmaxVals(2:end) - JmaxVals(1:end-1);

end