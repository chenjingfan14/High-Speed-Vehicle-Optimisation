function GlobalBest = gbestminroulette(nonDomFit,parsDominated,nonDomPos)
%% Roulette wheel global best selection
% Picks global best with respect to Pareto Front particles who dominate
% least number of particles. Therefore less particles dominated = higher
% chance of being selected (ie. outlier Pareto Front particles)

[dim,~] = size(nonDomFit);

num = (1:dim)';

dsort = sortrows([num,1./parsDominated],2);

partialSum = zeros(dim+1,1);
for i=2:dim+1
    partialSum(i) = dsort(i-1,2) + partialSum(i-1);
end

total = partialSum(end);

rtotal = rand*total;
ind = rtotal >= partialSum(1:end-1) & rtotal <= partialSum(2:end);

GlobalBest = nonDomPos(ind,:);

[check,~] = size(GlobalBest);

if check > 1
    ind = randi(check,1);
    GlobalBest = GlobalBest(ind,:);
end