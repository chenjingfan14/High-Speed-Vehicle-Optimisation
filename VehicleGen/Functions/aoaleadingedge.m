% function alteredwing = aoaleadingedge(del)

[row,col] = size(del);

for i = 1:floor(col/2)
    
    c1 = i;
    c2 = col - i + 1;
    
    [~,a] = max(del(:,c1));
    [~,b] = max(del(:,c2));
    
end