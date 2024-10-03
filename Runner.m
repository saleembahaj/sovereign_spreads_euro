%% ============================= Description  =============================

% This code produces the results from the VAR in "Sovereign Spreads in the 
% Euro Area: Cross Border Transmission and Macroeconomic Implications"

% The raw data is saved in VAR_database.xlsx. All the required functions
% and sub files are 

% The code produces Figures 4-7 from the paper.

% Written by Saleem Bahaj (Contact: saleembahaj@googlemail.com)
% Last updated: 28th December 2018

%% =====================LOAD DATA AND SPECIFICATION =======================

% This subfile loads in the data from excel the spreadsheet 
data_read;
% This subfile defines the specification
specification;

Iraw_d(isnan(Iraw_d(:,:,1)),:,:)=0; % Replace missing instrument observations with zero


%% ======================== Running estimation ============================      
[output] = PBVARX_HIERARCHICAL_FUN_COMP(Yraw_d,[],Iraw_d,varparaminputs);
%save VAR_output output varparaminputs Yraw_d Iraw_d -v7.3

%% =========================== GRAPH PRELIMS ==============================
% Labels
var_names=char('2 Year Yield','Industrial Production','Unemployment','Inflation','Financial Conditions','Fiscal Balance','German 2 Year Yield','Extra');
country_names2=char('Spain','Italy','Portugal','Ireland');

% Figure settings
colors = get(gca,'ColorOrder');
blue = colors(1,:);
red = colors(2,:);
add = 0.015;
add3 = 0.015;
add2 = add3/2;
add4 = add/4;
fstitle = 11;
fstic = 10;
fs = 10;

% Setup for IRFs
ihor=varparaminputs.ihor;

% Process output
for s=1:1
N=size(output.ac_draws,3);
gg=size(output.ac_draws,4);
M=size(output.ac_draws,2);
% Storage for the responses
Me=zeros(ihor,M,1);
Up=zeros(ihor,M,1);
Lo=zeros(ihor,M,1);
Up2=zeros(ihor,M,1);
Lo2=zeros(ihor,M,1);
Mec=zeros(ihor,M,1,N);
Upc=zeros(ihor,M,1,N);
Loc=zeros(ihor,M,1,N);
Upc2=zeros(ihor,M,1,N);
Loc2=zeros(ihor,M,1,N);

for i=1:gg
output.sir(:,:,:,i)=output.sir(:,:,:,i)/output.sir(1,1,1,i);
for k=1:M
if ismember(k,diffs)
output.sir(k,:,:,i)=cumsum(output.sir(k,:,:,i),3);  
end
end
for j=1:N
output.sirc(:,:,:,j,i)=output.sirc(:,:,:,j,i)/output.sirc(1,1,1,j,i);
for k=1:M
if ismember(k,diffs)
output.sirc(k,:,:,j,i)=cumsum(output.sirc(k,:,:,j,i),3);  
end
end
end
end
output.sir=sort(output.sir,4);
output.sirc=sort(output.sirc,5);
for j=1:M   
Me(:,j,1)=squeeze(output.sir(j,1,:,fix(gg*0.50)));        %median    
Up(:,j,1)=squeeze(output.sir(j,1,:,fix(gg*(1-irfp/2))));  %upper bound
Lo(:,j,1)=squeeze(output.sir(j,1,:,fix(gg*irfp/2)));      % lower bound
Up2(:,j,1)=squeeze(output.sir(j,1,:,fix(gg*(1-irfp2/2))));  %upper bound
Lo2(:,j,1)=squeeze(output.sir(j,1,:,fix(gg*irfp2/2)));      % lower bound
for i=1:N
Mec(:,j,1,i)=squeeze(output.sirc(j,1,:,i,fix(gg*0.50)));        %median    
Upc(:,j,1,i)=squeeze(output.sirc(j,1,:,i,fix(gg*(1-irfp/2))));  %upper bound
Loc(:,j,1,i)=squeeze(output.sirc(j,1,:,i,fix(gg*irfp/2)));      % lower bound
Upc2(:,j,1,i)=squeeze(output.sirc(j,1,:,i,fix(gg*(1-irfp2/2))));  %upper bound
Loc2(:,j,1,i)=squeeze(output.sirc(j,1,:,i,fix(gg*irfp2/2)));      % lower bound
end
end
end

%% ======================== Fig 4: Mean Country ===========================
impulse_select=1:1:M; % select the responses we are interested in.
imp_num=size(impulse_select,2);
t=1:ihor;            % Time periods for impulses
zline=zeros(size(t));
figure('Name',['Impulse Responses model type for mean coutry model'],'NumberTitle','off')
set(gcf,'Color',[1,1,1])
scale=Me(1,1);   
colors = get(gca,'ColorOrder');
blue = colors(1,:);
red = colors(2,:);
add = 0.015;
add3 = 0.015;
add2 = add3/2;
add4 = add/4;
fstitle = 11;
fstic = 10;
fs = 10;

for j=1:imp_num
set(gca,'FontSize',fstic,'fontname','times')
h = subplot(2,round(imp_num/2),j);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
jj=impulse_select(j);
[hl,hp]=boundedline(t,Me(:,jj),[Me(:,jj)/scale-Lo2(:,jj)/scale,Up2(:,jj)/scale-Me(:,jj)/scale],'-b',t,Me(:,jj)/scale,[Me(:,jj)/scale-Lo(:,jj)/scale,Up(:,jj)/scale-Me(:,jj)/scale],'-b','alpha',t,zline,0,'k');
axis tight;
title(var_names(jj,:),'fontname','times','FontSize',fstitle,'fontweight','normal');
end
xSize = 21; 
ySize = 12;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperPositionMode','auto')


%%  ================== Fig 5: Figure Multi Country IRF ====================
% Just countries
impulse_select=[1:1:6]; % select the responses we are interested in.
imp_num=size(impulse_select,2);
t=1:ihor;            % Time periods for impulses
zline=zeros(size(t));

for ss=1:1
figure('Name',['Country Specific Impulse Responses'],'NumberTitle','off')
set(gcf,'Color',[1,1,1])
for i=1:N
scale=Mec(1,1,1,i);    

for j=1:imp_num
h=subplot(N,imp_num,(i-1)*imp_num+j);
set(gca,'FontSize',fstic,'fontname','times')
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
jj=impulse_select(j);
[hl,hp]=boundedline(t,Mec(:,jj,1,i)/scale, [Mec(:,jj,1,i)/scale-Loc2(:,jj,1,i)/scale,Upc2(:,jj,1,i)/scale-Mec(:,jj,1,i)/scale],'-b',t,Mec(:,jj,1,i)/scale, [Mec(:,jj,1,i)/scale-Loc(:,jj,1,i)/scale,Upc(:,jj,1,i)/scale-Mec(:,jj,1,i)/scale],'-b','alpha',t,zline,0,'k');
axis tight;    

if i==1
title(var_names(jj,:),'fontname','times','FontSize',fstitle,'fontweight','normal');
end

if j==1
ylabel(country_names2(i,:),'fontname','times','FontSize',fstitle,'fontweight','normal');    
end 

end
end
end

%A4 is: 21.0 x 29.7cm
xSize = 29.7; 
ySize = 21;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperPositionMode','auto')


%% ================== Fig 6: Variance Decomposition =======================

Me_decomp=zeros(ihor,M,1);  
Up_decomp=zeros(ihor,M,1);
Lo_decomp=zeros(ihor,M,1);
Mec_decomp=zeros(ihor,M,1,N);
Upc_decomp=zeros(ihor,M,1,N);
Loc_decomp=zeros(ihor,M,1,N);
decomp_draws=sort(output.decomp_draws,3);
decompc_draws=sort(output.decompc_draws,4);

for j=1:M    
Me_decomp(:,j,1)=squeeze(decomp_draws(j,:,fix(gg*0.50)));        %median     
Up_decomp(:,j,1)=squeeze(decomp_draws(j,:,fix(gg*(1-irfp/2))));  %upper bound 
Lo_decomp(:,j,1)=squeeze(decomp_draws(j,:,fix(gg*irfp/2)));      % lower bound
for i=1:N
Mec_decomp(:,j,1,i)=squeeze(decompc_draws(j,:,i,fix(gg*0.50)));        %median     
Upc_decomp(:,j,1,i)=squeeze(decompc_draws(j,:,i,fix(gg*(1-irfp2/2))));  %upper bound 
Loc_decomp(:,j,1,i)=squeeze(decompc_draws(j,:,i,fix(gg*irfp2/2)));      % lower bound
end
end

impulse_select=[1:1:7];
for ss=1:1
figure('Name',['Variance decompositions for mean coutry model'],'NumberTitle','off')
set(gcf,'Color',[1,1,1])
for j=1:M
h=subplot(2,round(M/2),j);
set(gca,'FontSize',fstic,'fontname','times')
jj=impulse_select(j);
ph = get(h, 'pos');
ph(4) = ph(4) + add;
ph(3) = ph(3) + add2;
ph(2) = ph(2) - add4;
set(h, 'pos', ph);
boundedline(t,Me_decomp(:,jj),[Me_decomp(:,jj)-Lo_decomp(:,jj),Up_decomp(:,jj)-Me_decomp(:,jj)],'-b');
axis tight;
title(var_names(jj,:),'fontname','times','FontSize',fstitle,'fontweight','normal');  
end
end

%A4 is: 21.0 x 29.7cm

xSize = 21; 
ySize = 12;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperPositionMode','auto')





%% ===================== Figure 7: Historical Decomp ======================

sample=NaN(5*12+3,1);
for i=1:6
for j=1:12    
sample(j+12*(i-1),:)=datenum(2006+i, j, 1, 1, 1, 1);
end
end
sample(6*12+1,:)=datenum(2013, 1, 1, 1, 1, 1);
sample(6*12+2,:)=datenum(2013, 2, 1, 1, 1, 1);
sample(6*12+3,:)=datenum(2013, 3, 1, 1, 1, 1);

start=1;
[y_counter] = counter_fac(Yraw_d,output.betac_draws,output.bc_draws,output.ec_draws,varparaminputs.p,M,varparaminputs.detr,start);
y_counter=sort(y_counter,4);
med_counter=median(y_counter,4);
up_counter=zeros(size(med_counter));
lo_counter=zeros(size(med_counter));
N=size(output.ac_draws,3);
for i=1:N
up_counter(:,:,i)=squeeze(y_counter(:,:,i,fix(gg*(1-irfp/2))));
lo_counter(:,:,i)=squeeze(y_counter(:,:,i,fix(gg*(irfp/2))));
end

dt=1;

med_prem=squeeze(Yraw_d(start:end,dt,:)-med_counter(start:end,dt,:));
up_prem=squeeze(Yraw_d(start:end,dt,:)-up_counter(start:end,dt,:));
lo_prem=squeeze(Yraw_d(start:end,dt,:)-lo_counter(start:end,dt,:));
t=sample(start:(end));


% Graph of difference with error
figure('Name',['Counterfactual Analysis'],'NumberTitle','off')
set(gcf,'Color',[1,1,1])
for i=1:N
subplot(N/2,2,i)
set(gca,'FontSize',fstic,'fontname','times')
boundedline(t,med_prem(:,i),[med_prem(:,i)-lo_prem(:,i),up_prem(:,i)-med_prem(:,i)],'-b',t,zeros(size(t)),0,'-k');
datetick('x','mmm/yy','keepticks');
axis tight;
title(country_names2(i,:),'fontname','times','FontSize',fstitle,'fontweight','normal');
end    


xSize = 21; 
ySize = 12;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperPositionMode','auto')

% Graph vs level
figure('Name',['Counterfactual Analysis'],'NumberTitle','off')
set(gcf,'Color',[1,1,1])
for i=1:N
subplot(2,2,i),
set(gca,'FontSize',fstic,'fontname','times')
plot(t,Yraw_d(start:end,dt,i),'-r',t, med_counter(start:end,dt,i),'--b');
datetick('x','mmm-yy','keepticks'); 
axis tight;
title(country_names2(i,:),'fontname','times','FontSize',fstitle,'fontweight','normal');
if i==N
legend1=legend('Actual','Counterfactual');
set(legend1,'Location','Northwest','FontSize',9,'fontname','times');
legend boxoff
end
end

xSize = 21; 
ySize = 12;

set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
     ,'PaperPosition',[0 0 xSize ySize],'PaperPositionMode','auto')








