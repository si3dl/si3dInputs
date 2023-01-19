%% bathy4si3d.m
% This script converts point files defining lake bathymetry and shoreline
% and formats them to be read into si3D
% ------------------------------------
% based on code created by A. Cortes, 
% Last Modified - M. Swann Nov 2022
% -------------------------------------

% Notes - Input bathy and shoreline files must be specified in UTM
% coordinates

%% ... User inputs

function [f]=bathy4si3d_git(lakename,pathname,savepath,bathyfile,shoreline,dx)

% 
% lakename = 'Lago Llanq'; % lake name (needs to be 10 characters long)
% pathname = '/Users/micahswann/Dropbox/03-CHILELAGOS_TAHOE/32- Si3D/2_Llanquihue/2_Raw_Data/Bathy/1_csv'; % path to bathy files
% bathyfile = 'Llanquihue_Bathy_50m.csv'; % name of bathy point file
% shorelinefile = 'Llanquihue_Shoreline.csv'; % name of shoreline point file
% savepath = '/Users/micahswann/Dropbox/03-CHILELAGOS_TAHOE/32- Si3D/2_Llanquihue/3_Inputs/bathy'; % location to save files
% dx = 400; %dx = horizontal grid size for output file [m]
% 

%% ... Load raw data: Grid points and shoreline
cd(pathname);
bathy = load(bathyfile);

% Re-sizing the grid and adding the shore line
Z = bathy(:,1);
X = bathy(:,2);
Y = bathy(:,3);

% make depths negative
Z = -1*abs(Z);

% load shoreline points
shoreline = load(shoreline);% shoreline
Xs = shoreline(:,1);
Ys = shoreline(:,2);
Zs = zeros(length(Xs),1);

% plot raw data
f = figure('position',[306 105 1193 854]);
ax(1) = subplot(2,2,1);
plot(X,Y,'.')
hold on
plot(Xs,Ys,'r.')
title('raw data')

% Merge all vectors
Xa = [X;Xs];
Ya = [Y;Ys];
Za = [Z;Zs];
reference = unique(Zs)-0.001;

xmin = min(min(Xa)); xmax = max(max(Xa));
ymin = min(min(Ya)); ymax = max(max(Ya));

%% ... Regrid to dx resolution

[XX,YY] = meshgrid(xmin:dx:xmax,ymin:dx:ymax);
ZZ = griddata(Xa,Ya,Za,XX,YY);
[nx,ny] = size(ZZ);

% ... Add NaNs when elevation is equal to the reference
Zvector = ZZ(:);
dum = find(Zvector >= reference);
Zvector(dum) = nan;
ZZnew = reshape(Zvector,nx,ny);


%% ... Plotting
ax(2) = subplot(2,2,2);
contourf(XX, YY, ZZnew)
colorbar
hold on
plot(Xs,Ys,'r.')
title('contour plot')
% Depth
ZZd = ZZnew - reference;

ax(3) = subplot(2,2,3);
surf(XX, YY, ZZd)
shading interp
colorbar
%caxis([0 200])
hold on
plot(Xs,Ys,'k.');
title('surf plot')

%% redefining terms and obtain grid size
dy = dx; 
xg = XX; clear XX;
yg = YY; clear YY;
zg = flipud(ZZd);
dum = find(zg > -0.2); zg(dum) = -0.2; % define the minimum depth
dum = find(~isnan(zg));
zz = -99*ones(size(zg));
dum = find(~isnan(zg));
zz(dum)=zg(dum)*(-10); 
[ny,nx]=size(zz);

ax(4) = subplot(2,2,4);
p = pcolor(flipud(zg)); shading flat
chandle = colorbar;
title('final grid')

cd(savepath)
saveas(f,'bathy.fig')
%% create h file for si3dinput

fname =['h'];
fid = fopen(fname,'w+');

string1 = [lakename,'     (dx=  ',char(string(dx)),'m),   imx =  '];
l = length(string1);
if l > 37
    n = l-37;
    disp(['ERROR! Lake name is ',char(string(n)),' characters too long. Please change and rerun'])
    
else
     fprintf(fid,'%s',lakename,'      (dx=  ',char(string(dx)),'m),   imx =  ',num2str(nx),',jmx =  ',num2str(ny),',ncols = ',num2str(nx));   
    fprintf(fid,'\n');
    H1=['HV       V']; 
    for j=1:nx-1;
        H1=[H1,'   V'];
    end; 
    fprintf(fid,'%s\n',H1);
    fprintf(fid,'%s','     '); 
    for j=1:nx-1;
        fprintf(fid,' %4.0f',j+1); 
    end; 
    fprintf(fid,' %4.0f\n',nx+1);
    disp('****************************************************'); 
    disp('SPECIFICATIONS FOR NUMERICAL GRID');
    disp('****************************************************');
    disp(['SIZE  = ',num2str(size(zz))]);
    disp(' ****************************************************');
    dum = '%5d'; 
    for j=1:nx;
        dum=[dum,' %4.0f'];
    end
    dum=[dum,' \n'];
    for j=1:ny
        fprintf(fid,dum,ny-j+2, zz(j,:));
    end
    fclose(fid)

end
end