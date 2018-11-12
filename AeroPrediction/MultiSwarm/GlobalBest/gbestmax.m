function GlobalBest = gbestmax(parsDominated,nonDomCrowd,nonDomPos)
%% Global best selection
% Most dominating particle = global best

% Selects global best based on 1) no. particles dominated > 2) furthest
% away 1) from other non-dominated particles > 3) random of 2)
con = parsDominated == max(parsDominated);
ind = num(con);
if length(ind) > 1
    con = nonDomCrowd(ind) == max(nonDomCrowd(ind));
    if sum(con) > 1
        ind = datasample(ind(con > 0),1);
    else
        ind = ind(con);
    end
end

GlobalBest = nonDomPos(ind,:);