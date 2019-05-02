function [Kd, Edh, Ebn] = Engerer2Separation(Egh, latitudes, longitudes, time, averaging_period)
% The Engerer Separation Model was first described and validated in the
% paper: Engerer, N. A. 2015. Minute resolution estimates of the diffuse
% fraction of global irradiance for southeastern Australia. Solar Energy.
% 116, 215-237. Section 4. Development of new models.
% The paper presents 3 variations of the Engerer separation model: 1, 2 & 3
% The Engerer1 is fit to Non-Cloud Enhancement data, the Engerer2 to Cloud
% Enhancement data and Engerer3 to clear-sky data only.

% The predictor variables to fit B0-5 are based on:
% kt = clearness index
% zen = zenith angle
% ast = apparent solar time
% dktc = ktc - kt (ktc is the clear-sky clearness index = ghi,cs/e_exth.
% k_de = Proportion of kd attributable to cloud enhancement

% The original formulation was found to be not general enough to other
% parts of the world and at different temporal resolutions, and so in a new
% paper, the parameters have been updated. Bright, Jamie M. & Engerer,
% Nicholas A. 2019. Engerer2: Global re-parameterisation, update and
% validation of an irradiance separation model at different temporal
% resolutions. Journal of Renewable & Sustainable Energy. 11(2), xxx.

% The new parameters are defined in the paper as:
% ===================================================================================================================
% Parameter, Predictor,   Orig 1-min, New 1-min,  5-min       10-min      15-min      30-min,     1-hr,       1-day
% C,          - &         0.042336,   0.10562,    0.093936,   0.079965,   0.065972,   0.032675,   -0.0097539, 0.32726
% B0,         - &         -3.7912,    -4.1332,    -4.5771,    -4.8539,    -4.7211,    -4.8681,    -5.3169,    -9.4391
% B1,         K_t,         7.5479,    8.2578,     8.4641,     8.4764,     8.3294,     8.1867,     8.5084,     17.113
% B2,         AST,        -0.010036,  0.010087,   0.010012,   0.018849,   0.0095444,  0.015829,   0.013241,   0.13752
% B3,         zen,        0.003148,   0.00088801, 0.003975,   0.0051497,  0.0053493,  0.0059922,  0.0074356,  -0.024099
% B4,         dKtc,       -5.3146,    -4.9302,    -4.3921,    -4.1457,    -4.169,     -4.0304,    -3.0329,    6.6257
% B5,         Kde,        1.7073,     0.44378,    0.39331,    0.37466,    0.39526,    0.47371,    0.56403,    0.31419
% ===================================================================================================================

% INPUT PARAMETERS
% ghi   = Global horizontal irradiance [W/m^2]
% time = in the matlab datevec format [yyyy, mm, dd, HH, MM, SS] [UTC]
% latitudes = either a latitude per measurement or one for all.
% longitude = either a longitude per measurement or one for all.
% averaging_period = x in minutes. Must be 1, 5, 10, 15, 30, 60 or 1440 

% OUTPUT PARAMETERS
% Kd = diffuse fraction = Edh/ghi [fraction]
% Edh = diffuse horizontal irradiance [W/m^2]
% Ebn = direct/beam normal irradiance [W/m^2]
% 
% % EXMAPLE USAGE
% time = [2018,05,01,08,01,00;...
%         2016,03,05,10,15,00;...
%         2015,09,10,13,17,00;...
%         2014,01,20,20,18,00];
% Egh = [300;400;1000;600];
% longitudes=[-60;-5;20;170];
% latitudes=[-10;10;40;-40];
% averaging_period = 1;
% [Kd, Edh, Ebn] = Engerer2Separation(Egh, latitudes, longitudes, time, averaging_period)
% 


%% Input saftey checks
% check length of all data is equal
if length(unique([length(Egh),size(time,1)]))~=1
    error('Egh and time must be the same length of corresponding time steps')
end
% check class 
if ~isnumeric(Egh)
    error('Egh must be numeric')
end
if ~isnumeric(time)
    error('time must be numeric')
end

permissible_averaging_periods = [1, 5, 10, 15, 30, 60, 1440];
if sum(permissible_averaging_periods==averaging_period)~=1
    error('averaging_period must be numeric and be one of: 1, 5, 10, 15, 30, 60 or 1440')
end

%% Calculate solar angles and clear-sky irradiance 

% UTC timestep: 
% Decimal time (e.g 1.5 represents 01:30am)
decimal_time  = (time(:,4).*60 + time(:,5) + time(:,6)./60^2)./60; %where (:,4) is the hour, (:,5) is the minute, and (:,6) is the second.

% Day of year
doy = datenum( time(:,1:3) ) - datenum( [time(:,1) ones(size(time,1),2)] ) +1;  

% Equation of time
beta_eot = (360 / 365.242).*(doy-1); 
EoT = (0.258 .* cosd(beta_eot) - 7.416 .* sind(beta_eot) - 3.648 .* cosd(2 .* beta_eot) - 9.228 .* sind(2 .* beta_eot)); 

% Local solar noon
lsn = 12 - longitudes./15 -1.*(EoT./60) ;                                            

% Hour angle
hour_angle = (decimal_time-lsn).*15;  
hour_angle(hour_angle>=180)  = hour_angle(hour_angle>=180) - 360;
hour_angle(hour_angle<=-180) = 360 + hour_angle(hour_angle<=-180);
                                   
% Declination angle
phi_r = (2*pi/365) * (doy+(decimal_time/24)-1);
delta_r = 0.006918...
    - 0.399912*cos(phi_r) + 0.070257*sin(phi_r)...
    - 0.006758*cos(2.*phi_r) + 0.000907*sin(2.*phi_r)...
    - 0.002697*cos(3.*phi_r) + 0.001480*sin(3.*phi_r);

% Zenith angle
zen = acosd( sin(deg2rad(latitudes)) .* sin(delta_r) + cos(deg2rad(latitudes)) .* cos(delta_r) .* cosd(hour_angle) );

% Solar constant used by TJ model
Esc = 1366.1;

% Extraterrestrial irradiance
beta=(2.*pi.*doy)./365;
Eextn = Esc .* (1.00011 + 0.034221 .* cos(beta) + 0.00128 .* sin(beta) + 0.000719 .* cos(2.*beta) + 0.000077 .* sin(2.*beta) );

% Extraterrestrial horizontal irradiance, with night time correction.
E0h = Eextn.*cosd(zen);
E0h(zen>90)=0;

% The TJ clear-sky model
A = 1160 + 75 .* sin((360.*(doy-275))./365);
k = 0.174 + 0.035 .* sin((360.*(doy-100))./365);
C = 0.095 + 0.04 .* sin((360.*(doy-100))./365);
Ebncs = A .* exp(-k./cosd(zen));
Ebncs(zen>90) = 0;
Eghcs = Ebncs.*cosd(zen) + C.*Ebncs;
Eghcs(zen>90) = 0;

%% Calculate predictors.
% clearness index
Kt = Egh./E0h;

% clearness index of clearsky
Ktc = Eghcs./E0h;

% deviation of clearness from clear-sky-clearness
dKtc = Ktc - Kt;

% apparent solar time
AST =hour_angle./15 + 12;
% quality check asserted in the original Engerer2 R code.
AST(AST<0)=abs(AST(AST<0));

% Proportion of Kd attributable to cloud enhancement from Engerer2 R code.
cloud_enhancement_estimate = Egh - Eghcs;
cloud_enhancement_estimate(cloud_enhancement_estimate<0.015) = 0;
Kde = cloud_enhancement_estimate./Egh;

%% QC the predictors
% it is likely/possible for there to be inf values and NaN values.
Kt(isinf(Kt))=NaN;
Ktc(isinf(Ktc))=NaN;
dKtc(isinf(dKtc))=NaN;
AST(isinf(AST))=NaN;
Kde(isinf(Kde))=NaN;

Kt(zen>90)=NaN;
Ktc(zen>90)=NaN;
dKtc(zen>90)=NaN;
AST(zen>90)=NaN;
Kde(zen>90)=NaN;

%% Engerer version parametrisation
%     %original Engerer 2 variables.
%     % parameters
%     C = 4.2336E-2;
%     B0 = -3.7912;
%     B1 = 7.547948473;
%     B2 = -1.0035892E-2;
%     B3 = 3.147968E-3;
%     B4 = -5.31460674;
%     B5 = 1.707321551;    
    
% Define the parameters from the paper.
p.t1 = [0.105620,-4.13320,8.25780,0.0100870,0.000888010,-4.93020,0.443780];
p.t5 = [0.0939360,-4.57710,8.46410,0.0100120,0.00397500,-4.39210,0.393310];
p.t10 = [0.0799650,-4.85390,8.47640,0.0188490,0.00514970,-4.14570,0.374660];
p.t15 = [0.0659720,-4.72110,8.32940,0.00954440,0.00534930,-4.16900,0.395260];
p.t30 = [0.0326750,-4.86810,8.18670,0.0158290,0.00599220,-4.03040,0.473710];
p.t60 = [-0.00975390,-5.31690,8.50840,0.0132410,0.00743560,-3.03290,0.564030];
p.t1440 = [0.327260,-9.43910,17.1130,0.137520,-0.0240990,6.62570,0.314190];

parameters = p.(['t',num2str(averaging_period)]);
C = parameters(1);
B0 = parameters(2);
B1 = parameters(3);
B2 = parameters(4);
B3 = parameters(5);
B4 = parameters(6);
B5 = parameters(7);

% generalised logistic function.
Engerer2 = @(C, B0, B1, B2, B3, B4, B5, zen, AST, dKtc, Kde) ...
    C + (1-C) ./ (1 + exp(B0 + B1.*Kt + B2.*AST + B3.*zen + B4.*dKtc)) + B5.*Kde;

%% Calculate the outputs
% calculated the diffuse fraction
Kd = Engerer2(C, B0, B1, B2, B3, B4, B5, zen, AST, dKtc, Kde);

% quality control.
Kd(Kd>1)=1;
Kd(Kd<0)=0;

% calculate the diffuse irradiance
Edh = Egh.*Kd;

% calculate the direct normal irradiance
Ebn = (Egh - Edh)./cosd(zen);

%% Quality control
% this is particularly relevant at low zeniths.
Ebn(zen>90) = NaN;
Edh(zen>90) = NaN;
Kd(zen>90) = NaN;

end
