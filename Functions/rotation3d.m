% function rotation3d(theta,phi,psi)

points = 

if numel(theta) > 1
    
    theta = reshape(theta,1,1,[]);
    phi = reshape(phi,1,1,[]);
    psi = reshape(psi,1,1,[]);
end

sthe = sin(theta);
sphi = sin(phi);
spsi = sin(psi);

cthe = cos(theta);
cphi = cos(phi);
cpsi = cos(psi);

A = [cpsi .* cphi - cthe * sphi * spsi, cpsi * sphi + cthe * cphi * spsi, spsi * sthe;
    -spsi * cphi - cthe * sphi * cpsi, -spsi * sphi + cthe * cphi * cpsi, cpsi * sthe;
    sthe * sphi, -sthe * cphi, cthe];

