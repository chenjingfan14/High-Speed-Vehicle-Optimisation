function [L,H,LH] = fsiinterpolation(aero,struct)

aeroMesh = aero.Points;
aeroNodes = aero.centre;

aeroMesh = reshape(aeroMesh(:),[],3);
aeroNodes = reshape(aeroNodes(:),[],3);

con = all(isnan(aeroMesh),2);
aeroMesh(con,:) = [];
con = all(isnan(aeroNodes),2);
aeroNodes(con,:) = [];

[numelAeroMesh,~] = size(aeroMesh);
[numelAero,~] = size(aeroNodes);

oneAero = ones(numelAero,1);
oneAeroMesh = ones(numelAeroMesh,1);

%% Aero node/node and mesh/node coupling
for i = numelAeroMesh:-1:1
    
    for j = numelAero:-1:1
        
        vec = aeroMesh(i,:) - aeroNodes(j,:);
        aeroMeshRBF(i,j,:) = (vec(1)^2 + vec(2)^2 + vec(3)^2)^0.5;
    
    end
end

xyzAeroMesh = [oneAeroMesh, aeroMesh];

Ava = zeros(numelAeroMesh, numelAero + 4);
Ava(:,1:4) = xyzAeroMesh;
Ava(:,5:end) = aeroMeshRBF;

for i = numelAero:-1:1
    
    for j = numelAero:-1:1
        
        vec = aeroNodes(i,:) - aeroNodes(j,:);
        Ma(i,j,:) = (vec(1)^2 + vec(2)^2 + vec(3)^2)^0.5;
    
    end
end

PT = [oneAero, aeroNodes];
P = PT';
Minv = inv(Ma);

Caa = zeros(numelAero + 4);
Caa(1:4,5:end) = P;
Caa(5:end,1:4) = PT;
Caa(5:end,5:end) = Ma;

Mp = inv(P * Minv * PT);

aaPLY = Mp * P * Minv;
aaRBF = Minv - Minv * PT * Mp * P * Minv;

L = Ava * [aaPLY; aaRBF];

zero = zeros(size(L));

if nargin == 1
    
    L = [L zero zero; zero L zero; zero zero L];
    H = [];
    LH = [];
    return
end

%% Structural node/node and aero/struct nodes coupling

structNodes = struct.Nodes;
structNodes = reshape(structNodes(:),[],3);
[numelStruct,~] = size(structNodes);
oneStruct = ones(numelStruct,1);

for i = numelAero:-1:1
    
    for j = numelStruct:-1:1
        
        vec = aeroNodes(i,:) - structNodes(j,:);
        aeroStructRBF(i,j,:) = (vec(1)^2 + vec(2)^2 + vec(3)^2)^0.5;
    
    end
end

Aas = zeros(numelAero, numelStruct + 4);
Aas(:,1:4) = PT;
Aas(:,5:end) = aeroStructRBF;

for i = numelStruct:-1:1
    
    for j = numelStruct:-1:1
        
        % Euclidean distance between two nodes
        vec = structNodes(i,:) - structNodes(j,:);
        Ms(i,j,:) = (vec(1)^2 + vec(2)^2 + vec(3)^2)^0.5;
    
    end
end

PT = [oneStruct, structNodes];
P = PT';
Minv = inv(Ms);

Css = zeros(numelStruct + 4);
Css(1:4,5:end) = P;
Css(5:end,1:4) = PT;
Css(5:end,5:end) = Minv;

Mp = inv(P * Minv * PT);

asPLY = Mp * P * Minv;
asRBF = Minv - Minv * PT * Mp * P * Minv;

H = Aas * [asPLY; asRBF];

LH = L * H;

zero1 = zeros(size(H));
zero2 = zeros(size(LH));

L = [L zero zero; zero L zero; zero zero L];
H = [H zero1 zero1; zero1 H zero1; zero1 zero1 H];
LH = [LH zero2 zero2; zero2 LH zero2; zero2 zero2 LH];

end