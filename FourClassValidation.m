% This is a function that takes corresponding time series and performs
% an intense validation following the proceedure described by Chris A.
% Gueymard in his 2014 paper:
%
% Christian A. Gueymard. 2014. A review of validation methodologies and
% statistical performance indicators for modeled solar radiation data:
% Towards a better bankability of solar projects. Journal of Renewable and
% Sustainable Energy Reviews. Volume 39. Pages 1024-1034.
%
% -----------------------------------------------------------------------
%                               CONTENTS
% -----------------------------------------------------------------------
% The method is categorised into four classes labeled A, B C and D.
%
% A. Class A—indicators of dispersion
%  A.1. Mean bias difference (MBD)
%  A.2. Root mean square difference (RMSD).
%  A.3. Mean absolute difference (MAD)
%  A.4. Standard deviation of the residual (SD)
%  A.A. Coefficient of determination (R2)
%  A.6. Slope of best-fit line (SBF)
%  A.7. Uncertainty at 95% (U95)
%  A.8. t-statistic (TS)
%
% B Class B—indicators of overall performance.
%  B.1. Nash–Sutcliffe's efficiency (NSE)
%  B.2. Willmotts's index of agreement (WIA)
%  B.3. Legates's coefficient of efficiency (LCE)
%
% C. Class C—indicators of distribution similitude
%  C.1. Kolmogorov-Smirnov test integral (KSI)
%  C.2. OVER statistic of relative frequency of exeedence situations (OVER)
%  C.3 Combined Performance Index (CPI)
%
% D. Class D—visual indicators
%  D.1. Taylor diagram
%  D.2. Mutual information diagram
%  D.3. Boxplot Diagram
%
% -----------------------------------------------------------------------
%                               INPUTS
% -----------------------------------------------------------------------
% Observations and Predictions must be numerical and of the same size.
% Observations is treated as the correct value. They should be column
% vectors.
%
% class_selection must be a cell indicating which Class tests are to be
% performed. The default tests are all A, B, C and D. To indicate an
% individual test, set class_selection to {'B'} for example.
%
% -----------------------------------------------------------------------
%                               OUTPUTS
% -----------------------------------------------------------------------
% The output is a struct containing the results, depending on the selected
% class_selection. Within each struct contains the respective results.
% For example, validation_struct.MBD contains the mean bias difference of
% the validation.
%
% -----------------------------------------------------------------------
%                               EXAMPLES
% -----------------------------------------------------------------------
% sites = 500;
% time = (-4*pi():0.001:4*pi())';
% % all observations are the same for each site.
% Observations = 1000.*sin(time);
% Observations = repmat(Observations,[1,sites]);
% Predictions = zeros(size(Observations));
% % create random noise about the observations.
% for i = 1:sites
%     Predictions(:,i) = 1000.*sin(time) + 500.*rand().*cos(time) - 500.*rand().*sin(time).^2 + 500.*rand().*sin(time).^2;
% end
% % remove negative values and then the obs look sort of like irradiance
% Predictions(Predictions<0)=NaN;
% Observations(Observations<0)=NaN;
% 
% % validate using the four class system
% validation_struct = FourClassValidation(Observations,Predictions,{'A','B','C'});
% figure(1)
% plot(validation_struct.MBD)
% hold on
% plot(validation_struct.RMSD)
% plot(validation_struct.MAD)
% plot(validation_struct.SD)
% hold off
% legend('MBD','RMSD','MAD','SD')
% figure(2)
% plot(validation_struct.R2)
% hold on
% plot(validation_struct.SBF)
% plot(validation_struct.LCE)
% plot(validation_struct.NSE)
% hold off
% legend('R2','SBF','LSE','NSE')

function validation_struct = FourClassValidation(Observations,Predictions,class_selection)
%% preliminary checks and preprocessing
% default class selection to all
if ~exist('class_selection','var')
    class_selection = {'A','B','C','D'};
end

% checks
if (~isnumeric(Observations) || ~isnumeric(Predictions))
    error('Observations and Predictions must be of class numeric')
end
if size(Observations)~=size(Predictions)
    error('Observations and Predictions must be of equal size')
end
permissible_classes = {'A','B','C','D'};
for i = 1:length(class_selection)
    if max(strcmp(permissible_classes,class_selection{i}))==0
        error(['User defined class ''',class_selection{i},''' is not a valid class'])
    end
end

% pre-allocate memory for output
for i = 1:length(class_selection)
    switch class_selection{i}
        case 'A'
            validation_struct.MBD=zeros(1,size(Observations,2)).*NaN;
            validation_struct.nRMSD=zeros(1,size(Observations,2)).*NaN;
            validation_struct.RMSD=zeros(1,size(Observations,2)).*NaN;
            validation_struct.MAD=zeros(1,size(Observations,2)).*NaN;
            validation_struct.SD=zeros(1,size(Observations,2)).*NaN;
            validation_struct.R2=zeros(1,size(Observations,2)).*NaN;
            validation_struct.SBF=zeros(1,size(Observations,2)).*NaN;
            validation_struct.U95=zeros(1,size(Observations,2)).*NaN;
            validation_struct.TS=zeros(1,size(Observations,2)).*NaN;
        case 'B'
%             validation_struct.NSE=zeros(1,size(Observations,2)).*NaN;
            validation_struct.WIA=zeros(1,size(Observations,2)).*NaN;
            validation_struct.LCE=zeros(1,size(Observations,2)).*NaN;
        case 'C'
            validation_struct.KSI=zeros(1,size(Observations,2)).*NaN;
            validation_struct.OVER=zeros(1,size(Observations,2)).*NaN;
            validation_struct.CPI=zeros(1,size(Observations,2)).*NaN;
        case 'D'
    end
end
validation_struct.Om=zeros(1,size(Observations,2)).*NaN;
validation_struct.Pm=zeros(1,size(Observations,2)).*NaN;
validation_struct.N=zeros(1,size(Observations,2)).*NaN;

% loop through each time series, assuming each column is a unique site to
% validate.

for i = 1:size(Predictions,2)
    
    % extract to simple variable names
    O = Observations(:,i);
    P = Predictions(:,i);
    
    % removal of NaN and inf values
    not_nan_inds = (~isnan(O) & ~isnan(P) & ~isinf(O) & ~isinf(P));
    O = O(not_nan_inds);
    P = P(not_nan_inds);
    
    % define common usages
    Om = mean(O);
    Pm = mean(P);
    N = length(O);
    
    validation_struct.Om(1,i) = Om;
    validation_struct.Pm(1,i) = Pm;
    validation_struct.N(1,i) = N;
    
    %% Class A - indicators of dispersion
    % These are the indicators that the majority of readers should be most
    % familiar with. They are all expressed here in percent (of Om) rather than
    % in absolute units (W/m2 for irradiances, or MJ/m2 or kWh/ m2 for
    % irradiations) because non-expert stakeholders can much more easily
    % understand percent results. In any case, stating the value of Om in all
    % validation results allows the experts to convert back the percent figures
    % into absolute units if they so desire. Formulas in this section are well
    % established and do not need further references.
    
    if max(strcmpi(class_selection,'A')==1)
        
        %-------------------------------------------------------------------------
        % A.1  Mean bias difference (MBD)
        %     validation_struct.MBD = 100./Om .* sum(P-O); % Gueymard's version
        validation_struct.MBD(1,i) = 100./Om .* mean(P-O);
        % note that in the paper, the /N is missing, however, this must be
        % incorrect as no mean would be calculated without it.
        %-------------------------------------------------------------------------
        % A.2 Root mean square difference (RMSD) and normalised RMSD
        validation_struct.nRMSD(1,i) = 100./Om .* sqrt(mean((P-O).^2));
        validation_struct.RMSD(1,i) = sqrt(mean((P-O).^2));
        %-------------------------------------------------------------------------
        % A.3 Mean absolute difference (MAD)
        %     validation_struct.MAD = 100./Om .* sum(abs(P-O)); % Gueymard's version
        validation_struct.MAD(1,i) = 100./Om .* mean(abs(P-O));
        % note that in the paper, the /N is missing, however, this must be
        % incorrect as no mean would be calculated without it.
        %-------------------------------------------------------------------------
        % A.4 Standard deviation of the residual (SD)
        validation_struct.SD(1,i) = 100./Om .* sqrt( sum(N.*(P-O).^2)  - sum((P-O).^2) ) ./ N;
        %-------------------------------------------------------------------------
        % A.5 Coefficient of determination (R2)
        %     validation_struct.R2 = (  sum((P-Pm) .* (O-Om)) ./ sum((P-Pm).^2 .* (O-Om).^2)  ).^2; % Gueymard's version
        % This gives excessively large values, perhaps an order of 100 out?
        validation_struct.R2(1,i) = 1 - sum((O-P).^2) ./ sum((O-Om).^2);
        % Replacement taken from wikipedia.org/wiki/Coefficient_of_determination
        %     R2 does not indicate whether:
        %     - the independent variables are a cause of the changes in the dependent variable;
        %     - omitted-variable bias exists;
        %     - the correct regression was used;
        %     - the most appropriate set of independent variables has been chosen;
        %     - there is collinearity present in the data on the explanatory variables;
        %     - the model might be improved by using transformed versions of the existing set of independent variables;
        %     - there are enough data points to make a solid conclusion.
        %-------------------------------------------------------------------------
        % A.6 Slope of best-fit line (SBF)
        validation_struct.SBF(1,i) = sum((P-Pm).*(O-Om)) ./ sum((O-Om).^2);
        %-------------------------------------------------------------------------
        % A.7 Uncertainty at 95%
        validation_struct.U95(1,i) = 1.96 .* (validation_struct.SD(1,i).^2 + validation_struct.RMSD(1,i).^2).^0.5;
        %-------------------------------------------------------------------------
        % A.8 t-statistic (TS)
        validation_struct.TS(1,i) = ( (N-1) .* validation_struct.MBD(1,i).^2 ./ (validation_struct.RMSD(1,i).^2-validation_struct.MBD(1,i).^2) ).^0.5;
        %-------------------------------------------------------------------------
    end
    
    
    %% Class B - Indicators of overall performance
    % These are indicators that are less common in the solar field than those
    % of Class A. They convey relatively similar information as those of Class
    % A, with the cosmetic advantage that a higher value indicates a better
    % model.
    
    if max(strcmpi(class_selection,'B')==1)
        %-------------------------------------------------------------------------
        % B.1 Nash-Sutcliffe's efficiency (NSE)
        %         validation_struct.NSE(1,i) = 1 - sum((P-O).^2) ./ sum((O-Om).^2);
        %  THE NSE is removed as it is identical to  R2
        %-------------------------------------------------------------------------
        % B.2 Willmotts's index of agreement (WIA)
        validation_struct.WIA(1,i) = 1 - sum(P-O).^2 ./ sum(abs(P-Om) + abs(O-Om)).^2;
        %-------------------------------------------------------------------------
        % B.3 Lagates's coefficient of efficiency (LCE)
        validation_struct.LCE(1,i) = 1 - sum(abs(P-O)) ./ sum(abs(O-Om));
        %-------------------------------------------------------------------------
        % LSE and NSE vary between 1 for perfect agreement and -inf for complete
        % disagreement, whereas WIE varies only between 1 and 0.
        
    end
    
    %% Class C
    % The goal is to compare one or more cumulative frequency distribution
    % of modeled data to that of a reference dataset. Can one or more single
    % number provide a measure of the similitude between two or more
    % distributions? Substantial progress in that direction resulted from an
    % initial study by Polo et al. [1],who proposed to use the Kolmogorov–
    % Smirnov test when comparing different cumulative distribution functions
    % (CDFs), because of its advantage of being non- parametric and valid for
    % any kind of CDF. Espinar et al. [2] developed the method further, now
    % referring to it as the Kolmogorov–Smirnov test Integral (KSI)
    
    if max(strcmpi(class_selection,'C')==1)
        %-------------------------------------------------------------------------
        % C.1 Kolmogorov-Smirnov test Integral (KSI)
        % irradiance must be binned into x by intervals of n
        xbins = linspace(nanmin(O),nanmax(O),15);
        xmin = min(xbins);
        xmax = max(xbins);
        Od = histc(O,xbins)./N;
        Pd = histc(P,xbins)./N;
        % absolute difference between the two normalised distributions
        Dn = abs(Od - Pd);
        % pure function of N obtained from [3], though simplified to constant phi(N) \approx 1.63.
        Dc = 1.63/N^0.5;
        Ac = Dc .* (xmax - xmin);
        fun1 = @(x) interp1(xbins,Dn,x);
        validation_struct.KSI(1,i) = 100/Ac .* integral(fun1,xmin,xmax);
        % KSI is 0 if the two distributions being compared can be considered
        % identical in a statistical sense.
        
        %-------------------------------------------------------------------------
        % C.2 Relative frequency of exeedence situations (OVER)
        % The OVER test is the same as the KSI test, but only for those
        % bins that exceed the limit defined by Dc. This is useful when the
        % normalised distribution of modelled data points in specific bins
        % exceeds the critical limit that would make it statistically
        % undistinguisable from the reference distribution.
        Dnc = Dn;
        Dnc(Dnc<Dc) = 0;
        fun2 = @(x) interp1(xbins,Dnc,x);
        validation_struct.OVER(1,i) = 100/Ac .* integral(fun2,xmin,xmax);
        % OVER is 0 if the normalised distribution always remains below Dc.
        % OVER can be null indicating that the distribution of the
        % predictions generally respect those of the predictions.
        
        %-------------------------------------------------------------------------
        % C.3 Combined Performance Index (CPI)
        % The interest of CPI is that it combines conventional information
        % about dispersion and biase (through RMSD) with information about
        % distribution likenesses (through KSI and OVER), whilst maintaining a
        % high degree of discrimination between the different models. This
        % feature is of course highly desireable when comparing differnet
        % models of similar performance. This is arguably the most significant
        % statistic to compare different model performance.
        
        % The CPI requires the RMSD of the two indices, should class_section A
        % not have been performed, the RMSD must also be calculated
        if max(strcmpi(class_selection,'A'))==0
            % A.2 Root mean square difference (RMSD)
            validation_struct.RMSD = 100./Om .* ( sum(P-O).^2 ./ N ).^0.5;
        end
        % all values must be in percentages.
        validation_struct.CPI(1,i) = (validation_struct.KSI(1,i) + validation_struct.OVER(1,i) + 2.*validation_struct.RMSD(1,i))./4;
        
        %-------------------------------------------------------------------------
    end
    
end

%% Class D
% This category is completely different from the three previous ones
% because the goal here is to obtain a visualization rather than summary
% statistics in the form of a few numbers
if max(strcmpi(class_selection,'D')==1)
    % The first recommended plots for class D is a Taylor diagram detailed
    % by KE Taylor [4] that combines RMSD, SD and R2 into a single polar
    % diagram. It is ideal for comparing the performance of many different
    % models.
    
    % The second suggestion is a Mutual Information Diagram, which is a
    % revision of the Taylor diagram proposed by [5].
    
    % The box plot is another decent variation for demonstrating
    % performance at different sites.
    
    
end

end

%% References
%
% [1] Polo J, Zarzalejo LF, Ramirez L, Espinar B. Iterative filtering of
% ground data for qualifying statistical models for solar irradiance
% estimation from satellite data. Sol. Energy 2006; 80:240–7
%
% [2] Espinar B, Ramirez L, Drews A, Beyer HG, Zarzalejo LF, Polo J, et
% al. Analysis of different comparison parameters applied to solar
% radiation data from satellite and German radiometric stations. Sol Energy
% 2009;83:118–25.
%
% [3] Marsaglia G, Tsang WW, Wang J. Evaluating Kolmogorov's Distribution.
% J Stat Softw 2003;8:1–4
%
% [4] Taylor KE. Summarizing multiple aspects of model performance in a
% single diagram. J Geophys Res 2001;106D:7183–92.
%
% [5] Correa CD, Lindstrom P. The mutual information diagram for
% uncertainty visualization. Int J Uncertain Quantif 2013;3:187–201.
