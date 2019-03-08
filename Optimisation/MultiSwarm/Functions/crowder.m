function crowdingDistance = crowder(nonDomFitNorm,maxf,minf,n)
%% Crowding distance calculator
% Computes how close Pareto Front particles are to other particles on the
% PF. Used to delete PF particles if number of particles within PF >
% assigned PFmax

[dim,~] = size(nonDomFitNorm);

num = (1:dim)';

dist = zeros(dim,1);

for i=1:n
    fsort = sortrows([num,nonDomFitNorm(:,i)],2);
    dist(fsort(1,1)) = inf; 
    dist(fsort(end,1)) = inf;
    for j=2:dim-1
        dist(fsort(j),1) = dist(fsort(j),1)+(fsort(j+1,2) - fsort(j-1,2))/(maxf(i)-minf(i));
    end
end

crowdingDistance = dist;