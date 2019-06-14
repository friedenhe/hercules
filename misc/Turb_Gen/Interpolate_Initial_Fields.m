% This script interpolates the coarse fields to the refined fields

clear all;
clc;

% Scaling factor for U
% since the initial fields is for Ubar~=15, you can multiply
% a Uscaling to the U field to get the desired velocity field 
Uscaling=1.0

% domain size
Lx=4*pi;
Ly=4*pi/3;
Lz=2.0;

% coase grid numbers
nx1=32;
ny1=32;
nz1=32;

% refined grid numbers
% you need to change these numbers accordingly
nx2=256;
ny2=192;
nz2=128;

% calculate x,y,z for the coarse grids
x_nx1=linspace(0.0,Lx,nx1);
y_ny1=linspace(0.0,Ly,ny1);
z_nz1=get_z(0.0,Lz,nz1);

% calculate x,y,z for the refined grids
x_nx2=linspace(0.0,Lx,nx2);
y_ny2=linspace(0.0,Ly,ny2);
z_nz2=get_z(0.0,Lz,nz2);

% generate 3d mesh
[x1,y1,z1]=meshgrid(x_nx1,y_ny1,z_nz1);
[x2,y2,z2]=meshgrid(x_nx2,y_ny2,z_nz2);

for i=1:5
    
    if i==1
        inpath='./coarse_u_32_cubic.dat';
        outpath='./init_u.dat';
    end
    if i==2
        inpath='./coarse_v_32_cubic.dat';
        outpath='./init_v.dat';
    end
    if i==3
        inpath='./coarse_w_32_cubic.dat';
        outpath='./init_w.dat';
    end
    if i==4
        inpath='./coarse_p_32_cubic.dat';
        outpath='./init_p.dat';
    end
    if i==5
        inpath='./coarse_t_32_cubic.dat';
        outpath='./init_t.dat';
    end
    
    if (exist(inpath, 'file') ~=2)
        disp([inpath ' does not exist. Quit!'])
        return;
    end
    
    % read the coase flow fields
    fid=fopen(inpath,'r');
    var=fread(fid,nx1*ny1*nz1,'float64');
    
    % reshape the 1d var to 3d
    var3D1=reshape(var,nx1,ny1,nz1);
    
    % interpolate the coarse field to the refine one
    var3D2=interp3(x1,y1,z1,var3D1,x2,y2,z2,'spline');
    
    % reshape the 3d var to 1d one for output
    var1D2=reshape(var3D2,nx2*ny2*nz2,1);
    
    % output the fields
    outid=fopen(outpath,'w');
    % we may need some scaling for U
    if i==1
        fwrite(outid,var1D2*Uscaling,'float64'); % note that HERCULES needs double precision
    else
        fwrite(outid,var1D2,'float64'); % note that HERCULES needs double precision
    end
    
    disp(['Interpolating' outpath ' Finished!'])
    
    clear inpath outpath fid var var3D1 var3D2 var1D2 outid;
    
end



