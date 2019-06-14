% this function calculates z for channel flows
function [z] = get_z(zmin,zmax,nz)

zmean=0.5*(zmin+zmax);

alpha=1.1; % this parameter control the slop in the boundary layer
portion=0.5; % the grid is symmetric
a=(1+nz)*portion; % this is a parameter

for i=1:nz
    z(i)=zmean+(zmin-zmean)/tanh(2*alpha)*tanh(2*alpha*(i-a)/(1-a));
end

end
