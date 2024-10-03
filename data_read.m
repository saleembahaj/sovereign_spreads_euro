%% ============================= Description  =============================

% This brings in the data, organises it as a three dimensional data set,
% Options to take logs, demean and detrend the data are available.

% Should only be run as part of Runner.m

% Written by Saleem Bahaj (Contact: saleembahaj@googlemail.com)
% Last updated: 28th December 2018

%% Reading the data

regions=['ES';'IT';'PT';'IE'];
N=size(regions,1);
Y_raw_all=[];
for i=1:N   
raw_data = xlsread('VAR_Database.xlsx', regions(i,:), 'B3:I89');
raw_data(raw_data==-99999)=NaN;
Y_raw_all=cat(3,Y_raw_all,raw_data);
end
M=size(Y_raw_all,2);


% The Data has the following structure:
%1) Two year German bond yield
%2) IP - index 2005=100
%3) Private Sector Cost of Finance
%4) Two year bond yield
%5) Unemployment
%6) Primary fiscal balance
%7) Core CPI
%8) Instrument Benchmark


%% Transforming Data

% Options are:
% 0: No transformation
% 1: Log transformation
% 2: Log de-mean and linear detrend
% 3: Log de-mean and quadratic de-trend
% 4: Simple Difference - I throw away 1 observations from every series.
% 5: Log change - I throw away 1 observations from every series.
% 6: Annual log Change - I throw away 12 observations from every series.

d_options=zeros(size(raw_data,2),1)'; 
d_options(2)=6; % IP in annual growth rates
d_options(7)=6; % Core HICP in annual growth rates

d_options=kron(ones(N,1),d_options);

if not(size(d_options,2)==M)
error('Incorrect number of data options specified');   
end

t=size(Y_raw_all,1); % Number of time series observations
const=ones(t,1);     % Constant
trend=1:t;           % Trend
qtrend=(1:t).^2;     % Quadratic Trend 

Y_all=Y_raw_all;



% Transform Data
for i=1:N
for j=1:M
    
% do nothing    
if d_options(i,j)==0        

% take logs 
elseif d_options(i,j)==1  
Y_all(:,j,i)=log(Y_raw_all(:,j,i))*100;


% take logs and linear detrended
elseif d_options(i,j)==2    
y=log(Y_raw_all(:,j,i))*100;
y2=y;
y2(isnan(y),:)=[];
XXX=[const,trend'];
XXX2=XXX;
XXX2(isnan(y),:)=[];
res=y2-XXX2*inv(XXX2'*XXX2)*XXX2'*y2;
Y_all(isnan(y)~=1,j,i)=res;

% take logs and quadraric detrended
elseif d_options(i,j)==3   
y=log(Y_raw_all(:,j,i))*100;
y2=y;
y2(isnan(y),:)=[];
XXX=[const,trend',qtrend'];
XXX2=XXX;
XXX2(isnan(y),:)=[];
res=y2-XXX2*inv(XXX2'*XXX2)*XXX2'*y2;
Y_all(isnan(y)~=1,j,i)=res;

% simple first difference
elseif d_options(i,j)==4 
y=(Y_raw_all(:,j,i));
y2=[NaN(1,1);y(1:end-1,:)];
Y_all(:,j,i)=(y-y2);

% log first difference
elseif d_options(i,j)==5 
y=log(Y_raw_all(:,j,i));
y2=[NaN(1,1);y(1:end-1,:)];
Y_all(:,j,i)=(y-y2)*100;   

% log yoy difference
elseif d_options(i,j)==6 
y=log(Y_raw_all(:,j,i));
y2=[NaN(12,1);y(1:end-12,:)];
Y_all(:,j,i)=(y-y2)*100;

end
end
end

%Throwing away initial observations due to differencing
if max(max(d_options,[],1),[],2)==4||max(max(d_options,[],1),[],2)==5
Y_all=Y_all(2:end,:,:);
elseif max(max(d_options,[],1),[],2)==6
Y_all=Y_all(13:end,:,:);
end


% Throwing away observations to choose start
s=0; % setting s=0 gives 2007 as a start date 
Y_all=Y_all((s+1):end,:,:);

%save Y_all Y_all

