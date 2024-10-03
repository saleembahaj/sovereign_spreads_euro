%% ============================= Description  =============================
% This defines the specification for the PVAR model, generates data labels
% and builds the dataset.


% Should only be run as part of Runner.m

% Written by Saleem Bahaj (Contact: saleembahaj@googlemail.com)
% Last updated: 28th December 2018

%% ========================== PRELIMINARIES ===============================
% -------------- Define specification of the VAR model --------------------

varparaminputs.detr = 1;                    % 0: no determinstic, 1: intercepts, 
                                            % 2: intercept+linear trend, 3:intercept+quad trend   
varparaminputs.panel=[1,0,0];               % Set to 1 if you wish a deterministic component to 
                                            % be cross-section specific: [constant,linear,quadratic]       
varparaminputs.p = 2;                       % Number of lags on dependent variables
varparaminputs.q = [1,varparaminputs.p ];   % Number of lags on exogenous variables (redundant)
varparaminputs.ihor = 48;                   % Horizon to compute impulse responses
irfp=0.32;                                  % Bayesian Error Bands (inner)
irfp2=0.10;                                 % Bayesian Error Bands (outer)
varparaminputs.vu=11;                       % degrees of freedom on t-shocks (identification)
varparaminputs.nburn=50000;                 % size of burn-in
varparaminputs.thin=100;                    % thinning factor
varparaminputs.nsave=1000;                 % number of stored draws   
varparaminputs.initial=1;                   % Initialisation routine
                                            % 1= using country-by-country LS, 0=using pooled LS for initial values 


% ------------------------- Define variables ------------------------------
data_selection=[4 2 5 7 3 6 1]; % baseline
instrument_selection=[8];
M=size(data_selection,2);
%This Gives a data order of:
% Eurepo IP Prices bond_yield equity 10yr  
Ydat=[];

for jj=1:size(data_selection,2);
Ydat=[Ydat, Y_all(:,data_selection(jj),:)];     
end

Idat=[];

for jj=1:size(instrument_selection,2);
Idat=[Idat, Y_all(:,instrument_selection(jj),:)];     
end




variablenames=['  Sov Spread  ';
               '    Output    ';
               ' Unemployment ';
               '   Inflation  ';
               ' Cost of Fin. ';
               ' Fiscal Bal.  ';
               '  2 year OIS  ';
               '  Extra       ';
               '  10 year     '];

diffs=[];
%% ===================== Conditions for estimation ========================          
ihor=varparaminputs.ihor;


%------------------------------- Defining Panels --------------------------
% Selecting countries to be included
country_names=['ES';'IT';'PT';'IE'];
deficit=[1,2,3,4];
Yraw_d=[];
Iraw_d=[];
for ii=1:size(deficit,2);   
Yraw_d=cat(3,Yraw_d,Ydat(:,:,deficit(ii)));    
Iraw_d=cat(3,Iraw_d,Idat(:,:,deficit(ii)));  
end


Xraw=[];
Iraw_d=Iraw_d(1:end,:,:);
Yraw_d=Yraw_d(1:end,:,:);


if (max(max(isnan(Yraw_d))))==1
error('data contains NaNs')    
end 

