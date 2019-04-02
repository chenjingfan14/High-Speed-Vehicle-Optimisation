function partStruct = velocitydef(partStruct,run)
%% 

Vinf = num2cell(run.U);
[Vx,Vy,Vz] = deal(Vinf{:});

for ii=1:length(partStruct)
    
    % Initialise normal coords
    unitNorm = partStruct{ii}.unitNorm;
    
    dims = size(unitNorm);
    
    nx = unitNorm(:,:,1);
    ny = unitNorm(:,:,2);
    nz = unitNorm(:,:,3);

    [unitTang,unitSurf] = deal(zeros(dims));
    
    Tx = ny.*Vz - nz.*Vy;
    Ty = nz.*Vx - nx.*Vz;
    Tz = nx.*Vy - ny.*Vx;
    
    Tmag = (Tx.^2 + Ty.^2 + Tz.^2).^0.5;
    
    Sx = Ty.*nz - Tz.*ny;
    Sy = Tz.*nx - Tx.*nz;
    Sz = Tx.*ny - Ty.*nx;
    
    Smag = (Sx.^2 + Sy.^2 + Sz.^2).^0.5;
    
    % Only becomes unit at end
    unitTang(:,:,1) = Tx;
    unitTang(:,:,2) = Ty;
    unitTang(:,:,3) = Tz;
    
    unitSurf(:,:,1) = Sx;
    unitSurf(:,:,2) = Sy;
    unitSurf(:,:,3) = Sz;
    
    partStruct{ii}.unitTang = unitTang./Tmag;
    partStruct{ii}.unitSurf = unitSurf./Smag;
    
end