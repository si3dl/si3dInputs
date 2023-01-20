
%% initcond4si3d.m
% ----------------------------
% Generate initial condition file for Si3D 
% --------------------------------------
% Original code written by A.Cortes, February 2020
% Last Modified by M. Swann Nov 2022
%%--------------------------------------

% Notes
% Input - NxM matrix of temperatures in with N in time and M in depth
%       - spacing_method: constant - consant dz through water column
%                              - exponenital exponential 
% Ouputs - si3dinit.txt file


% Example
% pathname = '/Users/micahswann/Documents/GitHub/si3dInputs/_matllibrary_/Example'; % path to init cond file files
% fname = 'Llanquihue_T.mat';
% savepath = '/Users/micahswann/Dropbox/03-CHILELAGOS_TAHOE/32- Si3D/2_Llanquihue/3_Inputs/init_cond'; % location to save files
% TempProf = 'uniform'; % variable or uniform
% start_temp = 14.5; % starting initial temp if TempProf is UNIFORM
% spacing_method = 'constant'; % 
% dz = 2; % vertical grid size for CONSTANT layering
% tprof_path = '/Users/micahswann/Desktop/profile_example.txt';

function [f]= initcond4si3d_git(LakeName,fpath,savepath,SimStartDate,spacing_method,TempProf,dz,dzmin,dzmax);
    %% User Inputs
    
        % load initial conditions temp profile from tab deliminted text files
        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.DataLines = [1, Inf];
        opts.Delimiter = "\t";
        opts.VariableNames = ["depth", "temp"];
        opts.VariableTypes = ["double", "double"];
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        tprof = readtable(fpath, opts);
        z_obs = tprof.depth;
        t_obs = tprof.temp;
        n = length(z_obs);
        max_z = z_obs(end);
    
            if strcmp(TempProf,'uniform');
                t_obs = mean(t_obs,'omitnan');
            end

            %... add extra column for surface (= 0.5m temperature) for 1d interpolation
            str_temp = zeros(1,n+1);str_temp(1,[2:n+1]) = t_obs; str_temp(1,1) = t_obs(1);str_temp =str_temp';
            str_depth = zeros(1,n+1);str_depth(1,[2:n+1])= z_obs; str_depth = str_depth';

        if strcmp(spacing_method, 'constant')
            %...interpolate temperature at dz intervals
            
            depthi = (0:dz:max_z)';
            tempi = interp1(str_depth,str_temp,depthi); 
            z = -(depthi)'; z = flip(z);
            t = tempi; t = flip(t);
            
            %figure;plot(t,z)
     
            % ...Constant layering
            gridDZ = dz.*ones(1,1000);
            % ... Depths to top of layers in Si3D
            gridZ  = cumsum(gridDZ);
  
                    % ... Max depth (zl in si3d_inp.txt)
        zl = max_z; 
        dum = find(gridZ>=zl); km = dum(1); km1 = km + 2;  
        zlevel = [ -100 -100 gridZ(1:km)]'; % add -100 for boundary conditions
        [zu,dumu] = unique(z); tempu = t(dumu);
        [zs,dums] = sort(zu,'descend'); ts = tempu(dums);
        elseif strcmp(spacing_method, 'expo')


            n = [50]; % number of cells between dz_min and dz_ma≈
            gf = (dzmax/dzmin)^(1/(n - 1));
            k = 1:n;
            gridDZ = dzmin*gf.^(k-1);
            gridZ = cumsum(gridDZ);
            ii = length(gridZ);
            dum = gridZ(end);
            while dum<max_z
                dum = gridZ(ii)+dzmax;
                gridZ(ii+1) = dum;
                gridDZ(ii+1) = dzmax;
                ii = ii+1;
            end
%             gridZ(ii+1) = max_z;
%             gridDZ(ii+1) = max_z-gridZ(ii);
        

            depthi = (gridZ)';
            tempi = interp1(str_depth,str_temp,depthi); 
            z = -(depthi)'; z = flip(z);
            t = tempi; t = flip(t);
        
        % ... Max depth (zl in si3d_inp.txt)
        zl = max_z; 
        dum = find(gridZ>=zl); km = dum(1); km1 = km + 2;  
        zlevel = [ -100 -100 gridZ(1:km)]'; % add -100 for boundary conditions
        [zu,dumu] = unique(z); tempu = t(dumu);
        [zs,dums] = sort(zu,'descend'); ts = tempu(dums);
        end
        %% Generate data at depths defined in zlevel
        zz =  zlevel(2:end); zz(1) = 0;  
        zi = -(zz(1:end-1)+zz(2:end))/2.;
    
        ti = interp1(zs, ts, zi,'nearest','extrap');
    
        % Plot results
        % ============
        figure; plot(t_obs, -z_obs,'.',ti, zi,'-'); 
        legend('Real','Interpolated')
    
        figure('name','Initial Conditions Input');
        plot(ti,zi,'LineWidth',1.5);
        ylabel('Depth [m]');
        set(gca,'FontSize',14);
        xlabel('Temp [{\circ}C]');
        grid on

    %% ... Write Temperature data to input file for si3d

            % obtain date from doy to write in si3d file
        Dt =SimStartDate';
        m = month(Dt,'shortname');
        dom = day(Dt);
        str1  = strcat('Simulations starting on',{' '},m,{' '},string(dom), ', 00 h -',{' '});

    cd(savepath)
    fid = fopen('si3d_init.txt','wt+');
    fprintf(fid, '%s\n', 'Initial condition file for Si3D model  - ');
    fprintf(fid, '%s\n', 'Data from Lago Llanquihue database          - ');
    fprintf(fid, '%s\n', str1);
    fprintf(fid, '%s\n', 'Depths (m)   Temp (oC)                 - ');    
    fprintf(fid, '%s\n', 'Source: N/A           - ');
    fprintf(fid, '%s\n', '---------------------------------------- ');
    sixzeros = [ 0 0 0 0 0 0 ];
    fprintf(fid,'%10.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',zi(1), ti(1), sixzeros  );
    for j = 1:length(ti)
       fprintf(fid,'%10.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',zi(j), ti(j), sixzeros);
    end
    fprintf(fid,'%10.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',zi(j), ti(j), sixzeros);
    fclose(fid);
    
%% Write si3d_layer file

% ... Depths to top of layers in Si3D
z2 = (-zi);
for i = 2:length(z2);
    layertop(1,i) = (z2(i)+z2(i-1))./2;
end
% ... Max depth (zl in si3d_inp.txt)

dum = find(gridZ>=zl); km = dum(1); km1 = km + 2;  
zlevel = [ -100 -100 layertop]';
km1 = length(zlevel);
zl = max(z2);
% ... Plot
% figure; 
% subplot(3,1,1);plot(zlevel );xlabel('Layer');ylabel('Depth    (m)');
% hold on; plot([1 km1],[zl zl],'k', [km1 km1], [0 zl],'k');
% subplot(3,1,2);plot(zlevel );xlabel('Layer');ylabel('Thickness(m)');
% hold on; plot([1 km1],[gridDZ(km1) gridDZ(km1)],'k', [km1 km1], [0 gridDZ(km1)],'k');
% subplot(3,1,3);plot(gridDZ(1:km1-2), gridZ(1:km1-2));ylabel('Depth');xlabel('Thickness(m)');


%... Write zlevel data to input file for si3d
fid = fopen('si3d_layer.txt','wt+');
fprintf(fid, '%s\n', ['Depths to top of layers in Si3D Grid for Lake ',LakeName]);
fprintf(fid, '%s\n', '** used if ibathyf in si3d_inp.txt is set to < 0       ');
fprintf(fid, '%s\n', '------------------------------------------------------ ');
fprintf(fid,'%s \n', ['   km1   =        ', num2str(km1)]);
for j = 1:km1
   fprintf(fid,'%10.2f %10.4f \n',j, zlevel(j));
end
fclose(fid);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% dt2doy %%%
    function [dayofyear]  = dt2doy(dt)
    % converts datetime or datenum to doy of year 
        if isa(dt,'double') == 1
            dt = datetime(dt,'ConvertFrom','datenum');
        end
        yr = year(dt);
        dayofyear = zeros(length(dt),1);
        for i = 1:length(dt)
            str = datetime(yr(i),01,01)-1;
            str.TimeZone = dt.TimeZone;
            dur = dt(i)-str;
            dayofyear(i,1) = days(dur);
            clear dur
        end
    end

end