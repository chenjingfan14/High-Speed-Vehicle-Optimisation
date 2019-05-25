function [L,H,LH] = fsiinterpolation(RBFfun,aero,struct)

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
    
    diffMat = aeroMesh(i,:) - aeroNodes;
    X = (diffMat(:,1).^2 + diffMat(:,2).^2 + diffMat(:,3).^2).^0.5;
    aeroMeshRBF(i,:) = RBFfun(X);
end

aeroMeshRBF(isnan(aeroMeshRBF)) = 0;

xyzAeroMesh = [oneAeroMesh, aeroMesh];

Ava = zeros(numelAeroMesh, numelAero + 4);
Ava(:,1:4) = xyzAeroMesh;
Ava(:,5:end) = aeroMeshRBF;

for i = numelAero:-1:1
    
    diffMat = aeroNodes(i,:) - aeroNodes;
    X = (diffMat(:,1).^2 + diffMat(:,2).^2 + diffMat(:,3).^2).^0.5;
    Ma(i,:) = RBFfun(X);
end

Ma(isnan(Ma)) = 0;

PT = [oneAero, aeroNodes];
P = PT';
Minv = inv(Ma);

% Caa = zeros(numelAero + 4);
% Caa(1:4,5:end) = P;
% Caa(5:end,1:4) = PT;
% Caa(5:end,5:end) = Ma;

Mp = inv(P / Ma * PT);

aaPLY = Mp * P / Ma;

% diff = max(abs(test - aaPLY));

aaRBF = Minv - Minv * PT * Mp * P * Minv;

L = Ava * [aaPLY; aaRBF];

zero = zeros(size(L));

if nargin == 2
    
    L = [L zero zero; zero L zero; zero zero L];
    H = [];
    LH = [];
    
%     %% Check method
%     check = reshape(L * aeroNodes(:),[],3);
%     plotter(aero)
%     figure(gcf)
%     hold on
%     plot3(aeroNodes(:,1),aeroNodes(:,2),aeroNodes(:,3),'r.','DisplayName','Query Points');
%     plot3(check(:,1),check(:,2),check(:,3),'k.','DisplayName','Interpolated Points');
%     legend
%     hold off
    
    return
end

%% Structural node/node and aero/struct nodes coupling

structNodes = struct.Nodes;
structNodes = reshape(structNodes(:),[],3);
[numelStruct,~] = size(structNodes);
oneStruct = ones(numelStruct,1);

for i = numelAero:-1:1
    
    diffMat = aeroNodes(i,:) - structNodes;
    X = (diffMat(:,1).^2 + diffMat(:,2).^2 + diffMat(:,3).^2).^0.5;
    aeroStructRBF(i,:) = RBFfun(X);
end

aeroStructRBF(isnan(aeroStructRBF)) = 0;

%%
Aas = zeros(numelAero, numelStruct + 4);
Aas(:,1:4) = PT;
Aas(:,5:end) = aeroStructRBF;

for i = numelStruct:-1:1
    
    % Euclidean distance between two nodes
    diffMat = structNodes(i,:) - structNodes;
    X = (diffMat(:,1).^2 + diffMat(:,2).^2 + diffMat(:,3).^2).^0.5;
    Ms(i,:) = RBFfun(X);
end

Ms(isnan(Ms)) = 0;


PT = [oneStruct, structNodes];
P = PT';
Minv = inv(Ms);

% Css = zeros(numelStruct + 4);
% Css(1:4,5:end) = P;
% Css(5:end,1:4) = PT;
% Css(5:end,5:end) = Minv;

Mp = inv(P * Minv * PT);

asPLY = Mp * P * Minv;
asRBF = Minv - Minv * PT * Mp * P * Minv;

H = Aas * [asPLY; asRBF];

LH = L * H;

zero1 = zeros(size(H));
zero2 = zeros(size(LH));

% Aero node to aero corner
L = [L zero zero; zero L zero; zero zero L];
% Struct node to struct node
H = [H zero1 zero1; zero1 H zero1; zero1 zero1 H];
% Aero node to struct node
LH = [LH zero2 zero2; zero2 LH zero2; zero2 zero2 LH];

% %% Check method
% check = reshape(H * structNodes(:),[],3);
% plotter(aero)
% figure(gcf)
% hold on
% plot3(structNodes(:,1),structNodes(:,2),structNodes(:,3),'r.','DisplayName','Query Points');
% plot3(check(:,1),check(:,2),check(:,3),'k.','DisplayName','Interpolated Points');
% legend
% hold off
% 
% %% Check method
% check = reshape(H * aeroNodes(:),[],3);
% plotter(aero)
% figure(gcf)
% hold on
% plot3(aeroNodes(:,1),aeroNodes(:,2),aeroNodes(:,3),'r.','DisplayName','Query Points');
% plot3(check(:,1),check(:,2),check(:,3),'k.','DisplayName','Interpolated Points');
% legend
% hold off

end