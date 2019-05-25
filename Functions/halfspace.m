function val = halfspace(target,x,y)

% Halfspace search for lookup table data. find closest x to target and
% interpolate for y if required, otherwise return row indices

if size(target,2) ~= 1
    
    target = target(:);
end

dim = size(target);

nRows = length(x);

first = ones(dim);
last = first*nRows;
prevMiddle = zeros(dim);
stopCon = false(dim);
val = nan(dim);

while any(~stopCon)
    
    middle = floor((first + last)/2);
    
    lThan = target < x(middle);
    last(lThan) = middle(lThan) - 1;
    
    gThan = target > x(middle);
    first(gThan) = middle(gThan) + 1;
    
    con1 = prevMiddle == middle;
    
    if any(con1)
        
        stopCon(con1) = true;
    end
    
    prevMiddle = middle;
end

if nargin < 3
    
    val = middle;
else
    rows = middle;
    
    % Half space search can be +-1 away from actual closest lookup value
    % Therefore interpolate between rows,rows-1 and rows+1,rows and average
    
    con1 = rows == 1;
    con2 = rows == nRows;
    
    con = (con1 | con2);
    
    if any(con)
        
        [x0, x1, y0, y1] = deal(zeros(sum(con),1));
        
        r1 = rows(con1);
        r2 = rows(con2);
        
        x0(con1) = x(r1);
        x1(con1) = x(r1 + 1);
        
        y0(con1) = y(r1);
        y1(con1) = y(r1 + 1);
        
        x0(con2) = x(r2 - 1);
        x1(con2) = x(r2);
        
        y0(con2) = y(r2 - 1);
        y1(con2) = y(r2);
        
        val(con) = y0(con) + (target(con) - x0(con)).*((y1(con) - y0(con))./(x1(con) - x0(con)));
    end
    
    con = ~con;
    
    if any(con)
        
        r3 = rows(con);
        
        x0 = x(r3 - 1);
        x1 = x(r3);
        x2 = x(r3 + 1);
        
        y0 = y(r3 - 1);
        y1 = y(r3);
        y2 = y(r3 + 1);
        
        int1 = y0 + (target(con) - x0).*((y1 - y0)./(x1 - x0));
        int2 = y1 + (target(con) - x1).*((y2 - y1)./(x2 - x1));
        
        val(con) = (int1 + int2)/2;
    end
end