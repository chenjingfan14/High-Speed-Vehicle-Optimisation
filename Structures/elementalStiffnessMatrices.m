function k = elementalStiffnessMatrices(AE_L,E12_L3,E6_L2,GJ_L,E4_L,E2_L,Iy,Iz,lmn1,lmn2,lmn3)

[~,~,is3d] = size(AE_L);

if is3d == 1
    
    AE_L = permute(AE_L, [3 2 1]);
    E12_L3 = permute(E12_L3, [3 2 1]);
    E6_L2 = permute(E6_L2, [3 2 1]);
    GJ_L = permute(GJ_L, [3 2 1]);
    E4_L = permute(E4_L, [3 2 1]);
    E2_L = permute(E2_L, [3 2 1]);
    
    Iy = permute(Iy, [3 2 1]);
    Iz = permute(Iz, [3 2 1]);
    
end

lambdaMat = permute([lmn1, lmn2, lmn3],[2 3 1]);

l1 = lambdaMat(1,1,:);
l2 = lambdaMat(2,1,:);
l3 = lambdaMat(3,1,:);
m1 = lambdaMat(1,2,:);
m2 = lambdaMat(2,2,:);
m3 = lambdaMat(3,2,:);
n1 = lambdaMat(1,3,:);
n2 = lambdaMat(2,3,:);
n3 = lambdaMat(3,3,:);

l12 = l1.^2;
l22 = l2.^2;
l32 = l3.^2;
m12 = m1.^2;
m22 = m2.^2;
m32 = m3.^2;
n12 = n1.^2;
n22 = n2.^2;
n32 = n3.^2;

k1 = [(AE_L .* l12) + E12_L3.*(Iz .* l22 + Iy .* l32), (AE_L .* l1 .* m1) + E12_L3.*(Iz .* l2 .* m2 + Iy .* l3 .* m3), (AE_L .* l1 .* n1) + E12_L3.*(Iz .* l2 .* n2 + Iy .* l3 .* n3);...
    (AE_L .* l1 .* m1) + E12_L3.*(Iz .* l2 .* m2 + Iy .* l3 .* m3), (AE_L .* m12) + E12_L3.*(Iz .* m22 + Iy .* m32), (AE_L .* m1 .* n1) + E12_L3.*(Iz .* m2 .* n2 + Iy .* m3 .* n3);...
    (AE_L .* l1 .* n1) + E12_L3.*(Iz .* l2 .* n2 + Iy .* l3 .* n3), (AE_L .* m1 .* n1) + E12_L3.*(Iz .* m2 .* n2 + Iy .* m3 .* n3), (AE_L .* n12) + E12_L3.*(Iz .* n22 + Iy .* n32)];

k2 = [E6_L2.*(Iz - Iy) .* l2 .* l3, E6_L2.*(Iz .* l2 .* m3 - Iy .* l3 .* m2), E6_L2.*(Iz .* l2 .* n3 - Iy .* l3 .* n2);...
    E6_L2.*(Iz .* l3 .* m2 - Iy .* l2 .* m3), E6_L2.*(Iz - Iy) .* m2 .* m3, E6_L2.*(Iz .* m2 .* n3 - Iy .* m3 .* n2);...
    E6_L2.*(Iz .* l3 .* n2 - Iy .* l2 .* n3), E6_L2.*(Iz .* m3 .* n2 - Iy .* m2 .* n3), E6_L2.*(Iz - Iy) .* n2 .* n3];

k3 = [(GJ_L .* l12) + E4_L.*(Iy .* l22 + Iz .* l32), (GJ_L .* l1 .* m1) + E4_L.*(Iy .* l2 .* m2 + Iz .* l3 .* m3), (GJ_L .* l1 .* n1) + E4_L.*(Iy .* l2 .* n2 + Iz .* l3 .* n3);...
    (GJ_L .* l1 .* m1) + E4_L.*(Iy .* l2 .* m2 + Iz .* l3 .* m3), (GJ_L .* m12) + E4_L.*(Iy .* m22 + Iz .* m32), (GJ_L .* m1 .* n1) + E4_L.*(Iy .* m2 .* n2 + Iz .* m3 .* n3);...
    (GJ_L .* l1 .* n1) + E4_L.*(Iy .* l2 .* n2 + Iz .* l3 .* n3), (GJ_L .* m1 .* n1) + E4_L.*(Iy .* m2 .* n2 + Iz .* m3 .* n3), (GJ_L .* n12) + E4_L.*(Iy .* n22 + Iz .* n32)];

k4 = [(-GJ_L .* l12) + E2_L.*(Iy .* l22 + Iz .* l32), (-GJ_L .* l1 .* m1) + E2_L.*(Iy .* l2 .* m2 + Iz .* l3 .* m3), (-GJ_L .* l1 .* n1) + E2_L.*(Iy .* l2 .* n2 + Iz .* l3 .* n3);...
    (-GJ_L .* l1 .* m1) + E2_L.*(Iy .* l2 .* m2 + Iz .* l3 .* m3), (-GJ_L .* m12) + E2_L.*(Iy .* m22 + Iz .* m32), (-GJ_L .* m1 .* n1) + E2_L.*(Iy .* m2 .* n2 + Iz .* m3 .* n3);...
    (-GJ_L .* l1 .* n1) + E2_L.*(Iy .* l2 .* n2 + Iz .* l3 .* n3), (-GJ_L .* m1 .* n1) + E2_L.*(Iy .* m2 .* n2 + Iz .* m3 .* n3), (-GJ_L .* n12) + E2_L.*(Iy .* n22 + Iz .* n32)];

k2T = permute(k2,[2 1 3]);

k = [k1     k2     -k1      k2;...
    k2T    k3     -k2T     k4;...
    -k1    -k2      k1     -k2;...
    k2T    k4     -k2T     k3];