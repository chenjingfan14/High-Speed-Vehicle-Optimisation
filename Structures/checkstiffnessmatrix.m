function flag = checkstiffnessmatrix(K)

[~,~,dim] = size(K);

flag = zeros(dim,1);

for i = dim:-1:1
    
    k = K(:,:,i);
    
    notsym = ~issymmetric(k);
    % Seems to give very small < 0 values making it not positive definite
    [~,notposdef] = chol(k);
    
    % Second effort to remove these very small < 0 values
    if notposdef
        
        notposdef = any(eig(k) < -0.00001);
    end
    
    if notsym && notposdef
        
        flag(i) = 3;
        
    elseif notposdef
            
        flag(i) = 2;
        
    elseif notsym
        
        flag(i) = 1;
        
    end
end