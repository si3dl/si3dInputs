%% surfbc4si3d.m
% ----------------------------
% Generate surface boundary conditions file for si3d
% --------------------------------------
% Original code written by A.Cortes, February 2020
% Last Modified by M. Swann Jan 2023
%%--------------------------------------

% Notes
% Input - tab delinated text file specifying met variables in each column 
% in the following order:
%   Column 1 - dt - datetime in the format YYYY-MM-DD HH:mm:ss. Date should
%   be in local time for cloud cover calculations
%   Column 2 - attc [dimensionless] - light attenuation coefficient (1.7/secchi depth) can
%   be variable
%   Column 3 - Hsw [w m^-2] - incoming shorwtwave radiation 
%   Column 4 - Ta [C] - surface air temperature
%   Column 5 - Pa [kPa] - atmospheric pressure at lake surface
%   Column 6 - hr [%] - relative humidity as percent (e.g. 100 no 1 for   100%)
%   Column 7 - cw [dimensionless] - wind drag 9can be variable
%   Column 8 - WS [m/s] - wind speed at 10 m
%   Column 9 - WDir [degrees] - wind direction (0 degrees is north)

%% Surface Boundary condations (met data) for Beznar 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f]= surfbc4si3d(LakeName,metfile,savepath, surfbcType,SimStartDate,duration,timestep,lat,lon,lsm,hemisphere,spinup)
    reset(0)
    % import met data
    opts = delimitedTextImportOptions("NumVariables", 9);
    opts.DataLines = [2, Inf];
    opts.Delimiter = "\t";
    opts.VariableNames = ["dt", "attc", "Hsw", "Ta", "Pa", "hr", "cw", "WS", "WDir"];
    opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "dt", "InputFormat", "yyyy-MM-dd HH:mm:ss");
    tbl = readtable(metfile, opts);
    dt = tbl.dt;
    attc = tbl.attc;
    Hsw = tbl.Hsw;
    Ta = tbl.Ta;
    Pa = tbl.Pa*1000; % convert kPa to Pa
    hr = tbl.hr./100; % convert humidity to decimal form
    cw = tbl.cw;
    WS = tbl.WS;
    WDir = tbl.WDir;
    clear opts tbl

    % clip data to specified time frame
    dtS = SimStartDate;
    dtE = SimStartDate + duration/24/3600;
    ix = find(dt>=dtS & dt<=dtE);
    dt = dt(ix);
    attc = attc(ix);
    sw = Hsw(ix);
    Ta = Ta(ix);
    Pa = Pa(ix);
    rh = hr(ix);
    cw = cw(ix);
    ws = WS(ix);
    wd = WDir(ix);

    % get doy
    doyS = day(dtS,'dayofyear');
    doyE = day(dtE,'dayofyear');

    % get doy,day and hr
    yr = year(dt);uyr = unique(yr);
    dy = day(dt,'dayofyear');
    hr = hour(dt);
    mint = minute(dt);
    doy = dy+hr/24+mint/60/24;
    
    % plotting raw met data
    figure('name','Raw Data')
    ax(1) = subplot(6,1,1); plot(doy, Ta,'.'); set(gca,'xlim',[doyS doyE]);ylabel('Tair')
    ax(2) = subplot(6,1,2); plot(doy, rh,'.'); set(gca,'xlim',[doyS doyE]);set(gca,'ylim',[0 1]); ylabel('RH');
    ax(3) = subplot(6,1,3); plot(doy, ws,'.'); set(gca,'xlim',[doyS doyE]);ylabel('Wx Speed'); set(gca,'ylim',[0 15]);
    ax(4) = subplot(6,1,4); plot(doy, wd,'.'); set(gca,'xlim',[doyS doyE]);ylabel('WDir.'); set(gca,'ylim',[-2 360],'Ytick',[0:90:360])
    ax(5) = subplot(6,1,5); plot(doy, sw,'.'); set(gca,'xlim',[doyS doyE]); ylabel('SWin');
    ax(6) = subplot(6,1,6); plot(doy, Pa,'.'); set(gca,'xlim',[doyS doyE]); ylabel('Pa');
    linkaxes(ax,'x'); 
    xlabel('Day of year');

    %% ... Estimate cloudiness on a daily basis 

    % Calculate vapor pressure (mb) from Rh and Ta
    es = 6.11 * exp(17.3 * Ta ./ ( Ta + 237.3));
    ea = es .* rh; %%relative humidity as a fraction 
    
    % Astronomical calculations
    latitude = lat;                                 %Beznar
    Llm      = lon;                              % Local longitude (degrees)
    Lsm      = lsm ;                                % Local standard meridian (degrees)Sevilla(alcal� del rio)
    Hsc      = 1390;                                % Solar constant (W/m2)
    theta = latitude*pi/180;                        % latitude in radians
    I  = dy;                                       % Julian day (Jan.1 is julian day 1)
    r = 1 + 0.017 * cos((2*pi/365)*(186-I));        % Relative earth-sun distance
    d = 23.45 * pi/180 * cos((2*pi/365)*(172-I));   % Declination of the sun

    if strcmp(hemisphere,'S')
        d = -(d);                                       % flip declination if in southern hemisphere
    end
    dts = -1/15 * (Lsm-Llm);                        % Fraction of 15-degree increment Llm is east of Lsm
    value =         (sin(theta)*sin(d));
    value = value./ (cos(theta)*cos(d));
    tss = 12/pi * acos(-value) + dts + 12;          % Time of sunset
    tsu = -tss + 2 * dts + 24;                      % Time of sunrise
    gamma = zeros(size(hr)); 
    dum = find(hr>tsu & hr<tss); 
    gamma(dum) = 1;                                 % Correction factor
    dum1 = find(hr <=12 );                          % Hour angles - BEGINNING
    dum2 = find(hr > 12 );
    hb1  = pi/12*(hr-1-dts); 
    hb1(dum1) = hb1(dum1)+pi; 
    hb1(dum2) = hb1(dum2)-pi; 
    hb  = hb1;
    dum = find(hb1 > 2*pi); hb(dum) = hb(dum) + 2 * pi;
    dum = find(hb1 < 0   ); hb(dum) = hb(dum) - 2 * pi;
    dum1 = find(hr <=12 );                          % Hour angles - END
    dum2 = find(hr > 12 );
    he1  = pi/12*(hr-dts); 
    he1(dum1) = he1(dum1)+pi; 
    he1(dum2) = he1(dum2)-pi; 
    he  = he1;
    dum = find(he1 > 2*pi); he(dum) = he(dum) + 2*pi;
    dum = find(he1 < 0   ); he(dum) = he(dum) - 2*pi;
    Ho = Hsc./(r.^2).*(sin(theta).*sin(d)+12/pi* ...    % Extraterrestrial solar radiation
         cos(theta).*cos(d).*(sin(he)-sin(hb))).*gamma;
     
    % Assume an attenuation coefficient of approximately 0.76 (from p. 359
    % Martin & McCutcheon, 1999
    % ... Diminish Ho (external solar radiation by the attenuation coefficient)
    Ho = attc .* Ho; 
    dum = find(Ho<0.0); Ho(dum) = 0.0;
    
    % if missing Sw replace with Ho
    ix = find(isnan(sw));
    sw(ix) = Ho(ix);

    % Find out diurnal values of cloudiness
    cc = zeros(size(sw)); cc1 = cc;
    sro = sw;
    
    % MSwann - 5/4/20: srad sensor reading >0 at night which is
    % impacting cc calcs. Set all values <1 = 0
    idx = find(sro < 1); sro(idx) = 0;
    sra = Ho;

    cc  = real(sqrt((1-(sro./sra))./0.67)); % cloud cover
    io = isnan(cc); cc(io) = 0;
    cc1 = 1 - (sro ./ sra); % Fraction of solar radiation that reaches land
    
    % % ... Set limits to cloudiness
    dum = find(cc  > 1); cc (dum) = 1;
    dum = find(cc  < 0); cc (dum) = 0;
    dum = find(cc1 > 1); cc1(dum) = 1;
    dum = find(cc1 < 0); cc1(dum) = 0;
    
    % ... Select cloudiness
    Cl = cc;
    
    %% **** Estimation of ALBEDO ******************************* 
    % ... Use locally measured Qv and albedo 
    % a. First definition of albedo (Pivovarov, 1972)
    time = doy; 
    a0 = 0.03 + 0.01 * (0.5-cc).*(1.-sin(pi*(time-81)/183));
    declination = 0.4093 * sin( 2*pi*(time-79.75) / 365. );
    declination = -declination;
    sin_alpha = sin(lat*pi/180) * sin (declination) + ...
       cos(lat*pi/180) * cos (declination) .*	cos (pi / 12. * abs (hr - 12.));
    sin_alpha = max( 0.0, sin_alpha );
    alb = a0 ./ (a0 + sin_alpha);
    
    %[alb RefAng] = Albedo(doy,lat); % MacIntyre approach
    
    % a. Shortwave radiation
    Hswr   = sw .* (1 - alb); 
    SWin = sw;
    SWout = sw.*alb;
    figure('name','Shortwave Rad Comp');
    ax(1)= subplot(2,1,1);
    plot(doy,SWin);
    hold on
    plot(doy,SWout,doy,Hswr)
    legend('SWin','SWout','Hswr');
    ax(2) = subplot(2,1,2);
    plot(doy,alb)
    ylabel('Albedo (fraction)'); 
    hold on
    linkaxes(ax,'x');
    set(gca,'xlim',[doyS doyE])
    
    % figure comparing radiation
    clear ax
    figure('name','Measured and Calculated SWin Comparison');
    ax(1) = subplot(3,1,1); 
    plot(doy,sw,doy,Ho,'r'); set(gca,'xlim',[doyS doyE]);  
    legend('SWin-meas','SWin-cal');
    ylabel('Solar radiation (W/m^2)'); 
    ax(2) = subplot(3,1,2); 
    plot(doy,Cl,'o'); set(gca,'xlim',[doyS doyE]);  
    ylabel('Cloud cover (fraction)'); 
    ax(3) = subplot(3,1,3); 
    plot(doy,alb,'.-'); set(gca,'xlim',[doyS doyE]);
    ylabel('Albedo (fraction)'); 
    xlabel('Day of Year)');
    linkaxes(ax,'x'); 
    sgtitle('Lago Llanquihue')
    linkaxes(ax,'x')
    
    %rename data
    Ta1 = Ta;
    rh1   = rh; 
    Pa1   = Pa; 

    % calc u and v wind components
    u1 = -ws.*sin(wd*pi/180); 
    v1 = -ws.*cos(wd*pi/180); 

    % interpolate
    timestep = timestep/3600/24;
    doyR = [doyS:timestep:doyE];
    Ta   = interp1(doy,Ta1,doyR); 
    rh   = interp1(doy,rh1,doyR);
    u    = interp1(doy,u1 ,doyR);
    v    = interp1(doy,v1 ,doyR);
    cw   = interp1(doy,cw ,doyR);
    eta   = interp1(doy,attc ,doyR);
    cc   = interp1(doy ,Cl ,doyR);
    Hsw  = interp1(doy ,Hswr ,doyR); 
    Pa   = interp1(doy,Pa1,doyR);

    % Compute LWin
    LWin  = 0.937e-5 * 0.97 * 5.67e-8 * ((Ta + 273.16).^6) .* (1 + 0.17.*cc) ;

    % % .. Modify data during the specified spinup time
    dum = find(doyR < doyS+spinup/3600/24);
    u(dum) = 0.001;
    v(dum) = 0.001;

    
    figure('name','Si3D Input Data'); 
    ax(1) = subplot(7,1,1); plot(doyR, Ta); set(gca,'xlim',[doyS doyE]);ylabel('Tair')
    ax(2) = subplot(7,1,2); plot(doyR, rh); set(gca,'xlim',[doyS doyE]);set(gca,'ylim',[0 1]); ylabel('RH');
    ax(3) = subplot(7,1,3); plot(doyR, u); set(gca,'xlim',[doyS doyE]);ylabel('u');hold on; yline(0);
    ax(4) = subplot(7,1,4); plot(doyR, v); set(gca,'xlim',[doyS doyE]);ylabel('v');hold on; yline(0);
    ax(5) = subplot(7,1,5); plot(doyR, Hsw); set(gca,'xlim',[doyS doyE]); ylabel('SWnet');
    ax(6) = subplot(7,1,6); plot(doyR, LWin); set(gca,'xlim',[doyS doyE]); ylabel('LWin');
    ax(7) = subplot(7,1,7); plot(doyR, Pa); set(gca,'xlim',[doyS doyE]); ylabel('Pa');
    linkaxes(ax,'x'); 
    xlabel('Day of year (2021)');

%% ********************************** Write results for SI3D surfbc input file ********************************
    date = datetime(today,'ConvertFrom','datenum','Format','yyyy-MM-dd')'
    date = char(string(date));
    int = timestep*24*60; int = char(string(int));
    cd(savepath);
    r = length(doyR);
    fid=fopen('surfbc.txt','wt+');
    fprintf(fid,'%s \n', 'Surface boundary condition file for si3d model');
    fprintf(fid,'%s \n', [LakeName,' simulations']);
    fprintf(fid,'%s \n', ['Time is in given in hours from 00:00 hrs on julian day ', num2str(floor(doyS)), ',', num2str(year(1))]);                         
    fprintf(fid,'%s \n', '   Time in   // Data format is (10X,G11.2,...) Time attc Hsw Ta Pa hr Qlw cw ua va');    
    fprintf(fid,'%s \n', ['   ',int,'-min    // SOURCE = XXXXXXX']);
    fprintf(fid,'%s \n', [' intervals  (Note : file prepared on ', date,')']);
    fprintf(fid,'%s \n', ['   npts = ', num2str(r)]);
    value=doyR(1);
    for j=1:r
       a0= (doyR(j)-value)*24; % Time in hours
       a1= eta(j); 		% light attenuation coefficient
       a2= Hsw(j); 	    % Penetrative component of heat flux (albedo already taken into account)
       a3= Ta(j);  		% Air temperature
       a4= Pa(j); 	    % Atmospheric pressure
       a5= rh(j);		% relative humidty (fraction)

       if surfbcType == 2
          a6= cc(j); 		% cloud cover (fraction)
       elseif surfbcType == 3
          a6= LWin(j);     % LWin 
       end
       a7= cw(j) ; 		% **** Wind drag coefficient
       a8= u (j) ; 		% **** Wind speed in the EW direction
       a9= v (j) ; 		% **** Wind speed in the NS direction 
       if ( a4 >= 100000 )
            fprintf(fid,'%10.4f %10.4f %10.4f %10.4f %10.3f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
             a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
       else
           fprintf(fid,'%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
            a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
       end
    end
    fclose(fid);

end