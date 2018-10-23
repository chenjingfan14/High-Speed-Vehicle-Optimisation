function partStruct = velocitydef(partStruct,run)
%% 

Vinf = num2cell(run.U);
[Vx,Vy,Vz] = deal(Vinf{:});

for ii=1:length(partStruct)
    
    % Initialise normal coords
    unitNorm = partStruct(ii).unitNorm;
    
    [dim1,dim3] = size(unitNorm);
    
    X = 1:3:dim3;
    Y = X + 1;
    Z = Y + 1;
    
    nx = unitNorm(:,X);
    ny = unitNorm(:,Y);
    nz = unitNorm(:,Z);

    [unitTang,unitSurf] = deal(zeros(dim1,dim3));
    
    Tx = ny.*Vz - nz.*Vy;
    Ty = nz.*Vx - nx.*Vz;
    Tz = nx.*Vy - ny.*Vx;
    
    Tmag = (Tx.^2 + Ty.^2 + Tz.^2).^0.5;
    
    Sx = Ty.*nz - Tz.*ny;
    Sy = Tz.*nx - Tx.*nz;
    Sz = Tx.*ny - Ty.*nx;
    
    Smag = (Sx.^2 + Sy.^2 + Sz.^2).^0.5;
    
    unitTx = Tx./Tmag;
    unitTy = Ty./Tmag;
    unitTz = Tz./Tmag;
    
    unitSx = Sx./Smag;
    unitSy = Sy./Smag;
    unitSz = Sz./Smag;
    
    unitTang(:,X) = unitTx;
    unitTang(:,Y) = unitTy;
    unitTang(:,Z) = unitTz;
    
    unitSurf(:,X) = unitSx;
    unitSurf(:,Y) = unitSy;
    unitSurf(:,Z) = unitSz;
    
    partStruct(ii).unitTang = unitTang;
    partStruct(ii).unitSurf = unitSurf;
    
end