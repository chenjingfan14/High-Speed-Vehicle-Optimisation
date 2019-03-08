function theta = barycentric(p,x1,x2,x3,x4,y1,y2,y3,y4,t1,t2,t3,t4)

[~,dim2] = size(p);

x = p(1,:);
y = p(2,:);

W1(1,:) = ((y2 - y3).*(x - x3) + (x3 - x2).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
W2(1,:) = ((y3 - y1).*(x - x3) + (x1 - x3).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));

W1(2,:) = ((y4 - y1).*(x - x1) + (x1 - x4).*(y - y1))./((y4 - y1).*(x3 - x1) + (x1 - x4).*(y3 - y1));
W2(2,:) = ((y1 - y3).*(x - x1) + (x3 - x1).*(y - y1))./((y4 - y1).*(x3 - x1) + (x1 - x4).*(y3 - y1));

W3 = 1 - W1 - W2;

triangle = W1 >= 0 & W2 >= 0 & W3 >= 0;

con = all(triangle,1);

triangle(1,con) = false;

%% EXTRAPOLATE
col = 1:dim2;

if ~all(any(triangle,1))
    mat1 = [W1(1,:);W2(1,:);W3(1,:)];
    mat2 = [W1(2,:);W2(2,:);W3(2,:)];
    
    con1 = mat1 == inf;
    con2 = mat2 == inf;
    
    mat1(con1) = nan;
    mat2(con2) = nan;
    
    a = min(mat1,[],1);
    b = min(mat2,[],1);
    
    [choice,row] = max([a;b],[],1);
    
    bool = choice < 0;
    
    ID = row(bool) + (col(bool)-1)*2;
    
    triangle(ID) = true;
end

[~,row] = max(triangle,[],1);
ID = row + (col-1)*2;

%%

W1f = W1(triangle);
W2f = W2(triangle);
W3f = W3(triangle);

W1f = W1f';
W2f = W2f';
W3f = W3f';

theta(1,:) = W1f.*t1 + W2f.*t2 + W3f.*t3;
theta(2,:) = W1f.*t3 + W2f.*t4 + W3f.*t1;

theta = theta(ID);