% This function re-parameterises the Engerer separation model using the
% Engerer 2 (the best performing) as the initial state.
%
%               DEPENDENCIES
% FourClassValidation.m
%
%               INPUT PARAMETERS
% Egh   = Global horizontal irradiance [W/m^2]
% Eghcs = Clear-sky global horizontal irradiance [W/m^2]
% E0h   = Extraterrestrial horizontal irradiance [W/m^2]
% zen   = Zenith angle [degrees]
% ha    = Hour angle [degrees]
% version = string indicating which formulation of the Engerer Separation
%           model is desired. 'engerer1', 'engerer2', or 'engerer3'.
% plot_figures = true of false, whether or not to plot outputs.
% mode_choices = a cell where each entry is a different mode. A mode is a
%               string that produces an appropriate discretisation from
%               available data. Note that mode_choices use 'z' for zenith and 'k'
%               for clearness index!
%
%   % Example clear-sky index mode_choices
%   k_range = 0:0.1:1.1;
%   mode_choices = cell(length(k_range)-1,1);
%   for i=1:length(k_range)-1
%       mode_choices{i,1} = ['(k>=',num2str(k_range(i)),' & k<',num2str(k_range(i+1)),')'];
%   end
%
%   % Example zenith mode_choices
%   z_range = 0:10:90;
%   mode_choices = cell(length(z_range)-1,1);
%   for i=1:length(z_range)-1
%       mode_choices{i,1} = ['(z>=',num2str(z_range(i)),' & z<',num2str(z_range(i+1)),')'];
%   end
%
%                   OUTPUT PARAMETERS
% parameters = the new parameters for the Engerer separation
% errs = error metrics of new Kd to actual kd, struct.

%
% All the models were trained on 1-min data. This time averaging is likely
% crucial. All input data must therefore be of the same averaging
% procedure.
%
% For the parameterisation approach, we build four mode_choices. These are
% parameters for 4 scenarios defined as:
%
% parameters.clear = kc>=0.99 and zen<=70
% parameters.partiallyclear = kc<0.99 and kc>=0.7 and zen<=70
% parameters.cloudy = kc<0.7 and zen<=70
% parameters.lowsun = zen>70;
%
% % Training purposes:
% LoadData
% load('coeff_functions.mat');
% UpdateSolcastDataWithCorrections
% addpath('c:\userdata\brightj\Documents\Work\Projects\useful-matlab-tools\')
% Egh_solcast = satdata.himawari.solcast.ghi_new;
% Egh_obs = satdata.himawari.obs.GHIsum;
% Egh = Egh_solcast; % select which data to use!
% Kd = satdata.himawari.obs.DIFmeas./Egh_obs;
% % delete data where Kd==1.
% % inds = (Kd>0.995 & Kd<1.005);
% Eghcs = satdata.himawari.obs.ghics;
% zen = satdata.himawari.solcast.zenith;
% datevecs = satdata.himawari.solcast.time;
% Esc=1366.1;
% ndd = datenum(datevecs) - datenum([datevecs(:,1),ones(length(datevecs(:,1)),2),zeros(length(datevecs(:,1)),3)]);
% beta=(2.*pi.*ndd)./365;
% Eextn=Esc*(1.00011+0.034221*cos(beta)+0.00128*sin(beta)+0.000719*cos(2*beta)+0.000077*sin(2*beta));
% E0h = Eextn.*cosd(zen);
% E0h(zen>90)=0;
% [~, ~, ha] = latlon2solarzenazi_jb(satdata.himawari.lats, satdata.himawari.lons, satdata.himawari.solcast.time); %time = datevec in UTC
% version = 'engerer1';
% clearvars -except Kd Egh* Eghcs E0h zen ha version mode_choices
%
function [parameters, errors, original_errors, vars] = ParameteriseEngerer2Separation(Kd, Egh, latitudes, longitudes, time)
% Note that the inputs require the ACTUAL MEASURED Kd values.

%% Input saftey checks
% check length of all data is equal
if length(unique([length(Kd),length(Egh),size(time,1)]))~=1
    error('Egh, Kd and time must all be the same length of corresponding time steps')
end
% check class
if ~isnumeric(Kd)
    error('Kd must be numeric')
end
if ~isnumeric(Egh)
    error('Egh must be numeric')
end

% mode choices
% this can be used to discretise and reparameterise by different variables.
mode_choices = {'(~isnan(k))'};


%% Calculate solar angles and clear-sky irradiance 

% UTC timestep: 
% Decimal time (e.g 1.5 represents 01:30am)
decimal_time    = (time(:,4).*60 + time(:,5) + time(:,6)./60^2)./60; %where (:,4) is the hour, (:,5) is the minute, and (:,6) is the second.

% Day of year
doy = datenum( time(:,1:3) ) - datenum( [time(:,1) ones(size(time,1),2)] ) +1;  

% Equation of time
beta = (360 / 365.242).*(doy-1); 
EoT = (0.258 .* cosd(beta) - 7.416 .* sind(beta) - 3.648 .* cosd(2 .* beta) - 9.228 .* sind(2 .* beta)); 

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

% Solar constant
Esc=1366.1;
% doy = datenum(time) - datenum([time(:,1),ones(length(time(:,1)),2),zeros(length(time(:,1)),3)]);

% Extraterrestrial irradiance
beta=(2.*pi.*doy)./365;
Eextn= Esc .* (1.00011 + 0.034221 .* cos(beta) + 0.00128 .* sin(beta) + 0.000719 .* cos(2.*beta) + 0.000077 .* sin(2.*beta) );

% Extraterrestrial horizontal irradiane, with night time correction.
E0h = Eextn.*cosd(zen);
E0h(zen>90)=0;

% The TJ clear-sky model
A = 1160 + 75 .* sin((360.*(doy-275))./365);
k = 0.174 + 0.035 .* sin((360.*(doy-100))./365);
C = 0.095 + 0.04 .* sin((360.*(doy-100))./365);
Ebncs = A .* exp(-k./cosd(zen));
Eghcs = Ebncs.*cosd(zen) + C.*Ebncs;


%% Calculate predictors.

% clearness index
Kt = Egh./E0h;

% clearness index of clearsky
Ktc = Eghcs./E0h;

% clear-sky index
kc = Egh./Eghcs;

% deviation of clearness from clear-sky-clearness
dKtc = Ktc - Kt;

% apparent solar time
AST =hour_angle./15 + 12;

% Proportion of Kd attributable to cloud enhancement
Kde = nanmax(0,1-(Eghcs./Egh));

% Combine the training data

%% Initial parameters
% as the Engerer2 always seems to do better, we will always use the
% engerer2 paramterisation to initialise the estimates
C = 4.2336E-2;
B0 = -3.7912;
B1 = 7.5479;
B2 = -1.0036E-2;
B3 = 3.1480E-3;
B4 = -5.3146;
B5 = 1.7073;
x0 = [C, B0, B1, B2, B3, B4, B5];
dat = [Kt,AST,zen,dKtc,Kde];

%% Quality control
% We can only operate with actual values, no nan or inf.
idxs = (~isnan(sum(dat,2)) & ~isnan(Kd) & ~isinf(Kd) & ~isinf(sum(dat,2)) & Kt<1);

% specify the X dependent variables (e.g. all the inputs)
X = dat(idxs,:);

% extract zenith angle and clear-sky index for use in discretisation
z = zen(idxs);
k = Kt(idxs);

% sepcify the output corresponding Y values.
Y = Kd(idxs);

% %% Defined the Engerer2 formulation approach
% ComputeKd = @(b,X) b(1) + (1-b(1)) ./ (1 + exp(b(2) + b(3).*X(:,1) + b(4).*X(:,2) + b(5).*X(:,3) + b(6).*X(:,4))) + b(7).*X(:,5);
%
% %sum squared error cost function
% SSECF = @(b) sum((Y - ComputeKd(b,X)).^2);
%
% % find the new parameters with fminsearch
% [new_parameters, ~] = fminsearch(SSECF,x0);


%% Parameterise the mode_choices, USE z and k
% as this is a very complex optimisztion, we must relax the tolerances
% options = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','iter');
parameters_fminsearch = NaN(length(mode_choices),length(x0));
for m = 1:length(mode_choices)
    inds= eval(mode_choices{m});
    x = X(inds,:);
    y = Y(inds);
    
    % define the generalised logistic function
    Engerer2 = @(b,x) b(1) + (1-b(1)) ./ (1 + exp(b(2) + b(3).*x(:,1) + b(4).*x(:,2) + b(5).*x(:,3) + b(6).*x(:,4))) + b(7).*x(:,5);
    
    % define the error metric to minimise by
    SSE = @(b) sum((y - Engerer2(b,x)).^2);
    
    % derive new parameters.
    parameters_fminsearch(m,:) = fminsearch(SSE,x0);
end


%% Calculate the errors
close all
Kd_engerer2 = Engerer2(x0,X);
errs_engerer2 = FourClassValidation(Y,Kd_engerer2,{'A'});

Kd_fminsearch = NaN(size(Y));
for m = 1:length(mode_choices)
    inds = eval(mode_choices{m});
    Kd_fminsearch(inds) = Engerer2(parameters_fminsearch(m,:),X(inds,:));
end
errs_fminsearch = FourClassValidation(Y,Kd_fminsearch,{'A'});


%% Outputs
errors = errs_fminsearch;
original_errors = errs_engerer2;
parameters = parameters_fminsearch;
vars.Kt = Kt;

% %% Plot some results
% if exist('plot_figures','var')
%     % figure comparing the original Engerer2 performance to the new performance
%     f=figure('Name','compare to original','color','w');
%     f.Position(3)=f.Position(3)*2;
%     kd_vars = {'_engerer2','_fminsearch'};
%     
%     kd_titles= {'Original Engerer2','Optimised with fminsearch','Optimised with fsolve','Optimised with lsqnonlin'};
%     for i =1:length(kd_vars)
%         subplot(1,2,i)
%         hold on
%         % scatter(Y,Kd_engerer2,'.','MarkerEdgeAlpha',0.03)
%         dat = [eval(['Kd',kd_vars{i}]),Y];
%         bins = 120;
%         bin_range = linspace(0,1.2,bins);
%         [N_for_pcs] = hist3(dat,'Ctrs',{bin_range bin_range});
%         dat(dat>bin_range(end))=NaN;
%         [N,b] = hist3(dat,'Ctrs',{bin_range bin_range});
%         imagesc(b{1}([1 end]),b{2}([1 end]),100.*N./sum(sum(N)));
%         set(gca,'YDir','normal')
%         c2 = prctile(reshape(100.*N./sum(sum(N)),[numel(100.*N./sum(sum(N))),1]),99.9);
%         caxis([0 c2])
%         T = [255,255,255;255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0]./255; %http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=5
%         x = [0,linspace(0.00001,c2,size(T,1)-1)];
%         map = interp1(x./c2,T,linspace(0,1,bins));
%         colormap(gca,map)
%         plot(0:0.1:1.2,0:0.1:1.2,'k')
%         hold off
%         axis([0 1.2 0 1.2])
%         xlabel('Measured Kd')
%         ylabel('Estimated Kd')
%         errs = eval(['errs',kd_vars{i},';']);
%         text(0.05, 1.18, ['RMSE = ',num2str(errs.RMSD,3)])
%         text(0.05, 1.12, ['MBE = ',num2str(errs.MBD,3)])
%         text(0.05, 1.06, ['nRMSE = ',num2str(errs.nRMSD,3)])
%         text(0.05, 1.00, ['R2 = ',num2str(errs.R2,3)])
%         title(kd_titles{i})
%     end
%     %
%     
%     % plot each mode considered by the new expression
%     f=figure('Name','explore the mode_choices','color','w');
%     sps = numSubplots(length(mode_choices));
%     kd_titles= {'Original Engerer2','Optimised with fminsearch'};
%     for i = 1:length(mode_choices)
%         subplot(sps(1),sps(2),i)
%         eval(['x = Kd_fminsearch(',mode_choices{i},');']);
%         eval(['y = Y(',mode_choices{i},');']);
%         hold on
%         % scatter(x,y,'.','MarkerEdgeAlpha',0.08)
%         dat = [y,x];
%         bins = 120;
%         bin_range = linspace(0,1.2,bins);
%         [N_for_pcs] = hist3(dat,'Ctrs',{bin_range bin_range});
%         dat(dat>bin_range(end))=NaN;
%         [N,b] = hist3(dat,'Ctrs',{bin_range bin_range});
%         imagesc(b{1}([1 end]),b{2}([1 end]),100.*N./sum(sum(N)));
%         set(gca,'YDir','normal')
%         c2 = prctile(reshape(100.*N./sum(sum(N)),[numel(100.*N./sum(sum(N))),1]),99.9);
%         caxis([0 c2])
%         T = [255,255,255;255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0]./255; %http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=5
%         x = [0,linspace(0.00001,c2,size(T,1)-1)];
%         map = interp1(x./c2,T,linspace(0,1,bins));
%         colormap(gca,map)
%         plot(0:0.1:1.2,0:0.1:1.2,'k')
%         hold off
%         axis([0 1.2 0 1.2])
%         xlabel('Measured Kd')
%         ylabel('Engerer2 Estimated Kd')
%         text(0.05,1.1,['N = ',num2str(sum(eval(mode_choices{i})))])
%         title(mode_choices{i})
%     end
%     
%     
%     % plot the kd over kt
%     f=figure('Name','Look at Kt relationship','color','w');
%     f.Position=[16 45 1207 852];
%     var = {'Y','Kd_engerer2','Kd_fminsearch'};
%     kd_titles= {'Actual relationship','Original Engerer2','Optimised with fminsearch'};
%     sps = numSubplots(length(var));
%     for i = 1:length(var)
%         subplot(sps(1),sps(2),i)
%         dat = [eval(var{i}),k];
%         bins = 120;
%         bin_range = linspace(0,1.2,bins);
%         [N_for_pcs] = hist3(dat,'Ctrs',{bin_range bin_range});
%         dat(dat>bin_range(end))=NaN;
%         [N,b] = hist3(dat,'Ctrs',{bin_range bin_range});
%         hold on
%         imagesc(b{1}([1 end]),b{2}([1 end]),100.*N./sum(sum(N)));
%         set(gca,'YDir','normal')
%         c2 = prctile(reshape(100.*N./sum(sum(N)),[numel(100.*N./sum(sum(N))),1]),99.9);
%         caxis([0 c2])
%         T = [255,255,255;255,247,236;254,232,200;253,212,158;253,187,132;252,141,89;239,101,72;215,48,31;179,0,0;127,0,0]./255; %http://colorbrewer2.org/#type=sequential&scheme=BuPu&n=5
%         x = [0,linspace(0.000001,c2,size(T,1)-1)];
%         map = interp1(x./c2,T,linspace(0,1,bins));
%         colormap(gca,map)
%         hold off
%         axis([0 1 0 1.2])
%         xlabel('Kt')
%         ylabel('Kd estimate')
%         title(kd_titles{i})
%         text(0.05,1.1,['N = ',num2str(numel(eval(var{i})))])
%     end
%     
% end
end