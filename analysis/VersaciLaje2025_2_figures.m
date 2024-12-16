%% Figures from Versaci & Laje 2025
%
% RUN FIRST '../data/VersaciLaje2025_preprocess_data.m'
% set current directory to /analysis


clear;
close all;

% data paths (set current directory to /analysis)
data_dir = '../data/';
pow_dir = '../data/power analysis/';

% color map
NColors = 100;
indColors = [1 50 60 70];
col1 = morgenstemning(NColors);
col = col1(indColors,:);

fig_size_cm=[24 12];
fsize=12;
msize = 4;

titles={'NORMAL NoFBK','NORMAL WithFBK','HIGH NoFBK','HIGH WithFBK'};

confidence = 95;
alpha = 1-confidence/100;
probability = 1-alpha/2;


%% Figure 2: SDe vs CVF

% load data
load([data_dir 'CVF1INTRA.mat']);
load([data_dir 'STDASY.mat']);

Y=STDASY;
Y_text='SD_A (ms)';
X=CVF1INTRA;
X_text='CV_F (nondim.)';


h1=figure(1);


% plot data
for conditioni=1:4
    subplot(2,2,conditioni)
    plot(X(:,conditioni),Y(:,conditioni),'o','markersize',msize,'markerfacecolor',col(conditioni,:),'markeredgecolor',col(conditioni,:))
    hold on
    
    xlim([0.07 0.4]);
    ylim([10 31]);
    if conditioni==1
        set(gca,'xtick',[0 0.1 0.2 0.3 0.4],'xticklabel',[],'ytick',[10 15 20 25]);
    elseif conditioni==2
        set(gca,'xtick',[0 0.1 0.2 0.3 0.4],'xticklabel',[],'ytick',[10 15 20 25],'yticklabel',[]);
    elseif conditioni==3
        set(gca,'xtick',[0 0.1 0.2 0.3 0.4],'xticklabel',[0 0.1 0.2 0.3 0.4],'ytick',[10 15 20 25]);
    elseif conditioni==4
        set(gca,'xtick',[0 0.1 0.2 0.3 0.4],'xticklabel',[0 0.1 0.2 0.3 0.4],'ytick',[10 15 20 25],'yticklabel',[]);
    end
    set(gca,'fontsize',fsize);
    set(gca,'box','off') %Remove y right ticks
    title([titles{:,conditioni}],'fontsize',fsize,'fontweight','normal')               
end


% fit Linear Mixed Model (LMM)
Y=STDASY; X=CVF1INTRA;
YY=[Y(:,1);Y(:,2);Y(:,3);Y(:,4)];
XX=[X(:,1);X(:,2);X(:,3);X(:,4)];
CONf=[zeros(22,1);ones(22,1);zeros(22,1);ones(22,1)];
HIGH=[zeros(44,1);ones(44,1)];
subject=[[1:22]';[1:22]';[23:44]';[23:44]'];

tbl=table(subject,YY,CONf,HIGH,XX);
% lme = fitlme(tbl,'YY ~ 1 + CONf + HIGH + XX + CONf:HIGH + XX:CONf + XX:HIGH + (1|subject)');
lme = fitlme(tbl,'YY ~ 1 + CONf + HIGH + XX + CONf:HIGH + XX:CONf + XX:HIGH + (1|subject)','DummyVarCoding','effects');
lme_anova = anova(lme);


% slope and intercept for every condition
b(1)=lme.Coefficients.Estimate(1);
b(2)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2);
b(3)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(3);
b(4)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2)+lme.Coefficients.Estimate(3)+lme.Coefficients.Estimate(5);

m(1)=lme.Coefficients.Estimate(4);
m(2)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6);
m(3)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(7);
m(4)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6)+lme.Coefficients.Estimate(7);


% R_squared calculation
for i=1:4    
    y_mean=mean(Y(:,i));
    SS_tot=sum((Y(:,i)-y_mean).^2);
    SS_res=0;
    for j=1:length(Y(:,1))
        temp=(m(i)*X(j,i)+b(i)-Y(j,i))^2;
        SS_res=SS_res+temp;
    end
    RSquared_LMM=1-(SS_res/SS_tot);
    p=2; %number of parameters of the model (slope and intercept)
    L=length(Y(:,1)); %sample size
    RSquared_LMM_adj=1-(1-RSquared_LMM)*(L-1)/(L-p);
end


% confidence intervals of regression
% ref: https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals
alpha/2;
n=size(X,1);
df=n-2; %degrees of freedom
t_score=tinv(probability,df);

% plot curve fit and calculate and plot confidence intervals
for conditioni=1:4
    %LMM fit
    subplot(2,2,conditioni);
    
    slope=m(conditioni);
    intercept=b(conditioni);
    x=X(:,conditioni);
    y=Y(:,conditioni);
    xFit = linspace(min(x),max(x), 20);%300
    yFit = slope*xFit+intercept;
    yFit_onlyExperimentalValues=slope*x+intercept;
    hold on
    plot(xFit,yFit,'-','MarkerSize',1,'LineWidth',2,'Color','k')
    MSE=sum((y-yFit_onlyExperimentalValues).^2)/(n-2);
    SSx=sum((x-mean(x)).^2);
    DIF2=(xFit-mean(x)).^2;
    error=t_score*sqrt(MSE*(1/n+DIF2./SSx));
    upper=yFit+error;
    lower=yFit-error;
    hold on
    patch([xFit fliplr(xFit)], [lower fliplr(upper)], [0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none')
end


h1_ax = axes(h1,'visible','off'); 
h1_ax.Title.Visible = 'on';
h1_ax.XLabel.Visible = 'on';
h1_ax.YLabel.Visible = 'on';
xlabel(h1_ax,X_text);
ylabel(h1_ax,Y_text);
set(gca,'fontsize',fsize);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h1);


disp([newline '<strong>FIGURE 2: SDe vs CVF</strong>']);
disp(['R^2=' num2str(round(100*RSquared_LMM,4)/100)]);
disp(lme);
disp(lme_anova);
print(h1,'-dpdf','-bestfit','./figure2_SDA_vs_CVF.pdf');



% printing as .eps doesn't work well in my version of Matlab, so I have to print as .pdf and convert afterwards
% Prevent font rasterization when converting to eps:
% Open pdf file in Inkscape, identify and select the transparent objects and delete them, then save as eps.

% https://stackoverflow.com/questions/72981242/converting-pdf-to-eps-without-rasterizing-or-changing-fonts
% Find/Replace (Ctrl+F) to search objects with string "clipPath" or "path" etc and
% with 'Search option = Properties'. Then I open the Objects Tab (Menu Object->
% Objects...) and use that to delete each transparent object generated by Powerpoint.




%% Figure 4: SDF vs F

% load data
load([data_dir 'F1.mat']);
load([data_dir 'STDF1.mat']);

Y=STDF1;
Y_text='SD_F (a.u.)';
X=F1;
X_text='F (a.u.)';


h2=figure(2);


% plot data
for conditioni=1:4
    subplot(2,2,conditioni)
    plot(X(:,conditioni),Y(:,conditioni),'o','markersize',msize,'markerfacecolor',col(conditioni,:),'markeredgecolor',col(conditioni,:))
    hold on
    
    ylim([0 0.8]);
    xlim([0 6]);
    if conditioni==1
        set(gca,'xtick',[0 2 4],'xticklabel',[],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[0 "" 0.2 "" 0.4 "" 0.6]);
    elseif conditioni==2
        set(gca,'xtick',[0 2 4],'xticklabel',[],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[]);
    elseif conditioni==3
        set(gca,'xtick',[0 2 4 6],'xticklabel',[0 2 4 6],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[0 "" 0.2 "" 0.4 "" 0.6]);
    elseif conditioni==4
        set(gca,'xtick',[0 2 4 6],'xticklabel',[0 2 4 6],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[]);
    end
    set(gca,'fontsize',fsize);
    set(gca,'box','off') %Remove y right ticks    
    title([titles{:,conditioni}],'fontsize',fsize,'fontweight','normal')
end


% fit Linear Mixed Model (LMM)
Y=STDF1; X=F1;
YY=[Y(:,1);Y(:,2);Y(:,3);Y(:,4)];
XX=[X(:,1);X(:,2);X(:,3);X(:,4)];
CONf=[zeros(22,1);ones(22,1);zeros(22,1);ones(22,1)];
HIGH=[zeros(44,1);ones(44,1)];
subject=[[1:22]';[1:22]';[23:44]';[23:44]'];

tbl=table(subject,YY,CONf,HIGH,XX);
lme = fitlme(tbl,'YY ~ 1 + CONf + HIGH + XX + CONf:HIGH + XX:CONf +XX:HIGH+ (1|subject)');
lme_amova = anova(lme);

% slope and intercept for every condition
b(1)=lme.Coefficients.Estimate(1);
b(2)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2);
b(3)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(3);
b(4)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2)+lme.Coefficients.Estimate(3)+lme.Coefficients.Estimate(5);

m(1)=lme.Coefficients.Estimate(4);
m(2)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6);
m(3)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(7);
m(4)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6)+lme.Coefficients.Estimate(7);

% R_squared calculation
for i=1:4    
    y_mean=mean(Y(:,i));
    SS_tot=sum((Y(:,i)-y_mean).^2);
    SS_res=0;
    for j=1:length(Y(:,1))
        temp=(m(i)*X(j,i)+b(i)-Y(j,i))^2;
        SS_res=SS_res+temp;
    end
    RSquared_LMM=1-(SS_res/SS_tot);
    p=2; %number of parameters of the model (slope and intercept)
    L=length(Y(:,1)); %sample size
    RSquared_LMM_adj=1-(1-RSquared_LMM)*(L-1)/(L-p);
end


% confidence intervals of regression
% ref: https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals
n=size(X,1);
df=n-2; %degrees of freedom
t_score=tinv(probability,df);

% plot curve fit and calculate and plot confidence intervals
for conditioni=1:4
    %LMM fit
    subplot(2,2,conditioni);
    
    slope=m(conditioni);
    intercept=b(conditioni);
    x=X(:,conditioni);
    y=Y(:,conditioni);
    xFit = linspace(min(x),max(x), 20);
    yFit = slope*xFit+intercept;
    yFit_onlyExperimentalValues=slope*x+intercept;
    hold on
    plot(xFit,yFit,'-','MarkerSize',1,'LineWidth',2,'Color','k')
    MSE=sum((y-yFit_onlyExperimentalValues).^2)/(n-2);
    SSx=sum((x-mean(x)).^2);
    DIF2=(xFit-mean(x)).^2;
    error=t_score*sqrt(MSE*(1/n+DIF2./SSx));
    upper=yFit+error;
    lower=yFit-error;
    hold on
    patch([xFit fliplr(xFit)], [lower fliplr(upper)], [0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none')
end



h2_ax = axes(h2,'visible','off');
h2_ax.Title.Visible = 'on';
h2_ax.XLabel.Visible = 'on';
h2_ax.YLabel.Visible = 'on';
xlabel(h2_ax,X_text);
ylabel(h2_ax,Y_text);
set(gca,'fontsize',fsize);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h2);


disp([newline '<strong>FIGURE 3: SDF vs F</strong>']);
disp(['R^2=' num2str(round(100*RSquared_LMM,4)/100)]);
disp(lme);
disp(lme_anova);
print(h2,'-dpdf','-bestft','./figure4_SDF_vs_F.pdf');



%% Figure 5: CVF vs F

% load data
load([data_dir 'F1.mat']);
load([data_dir 'CVF1INTRA.mat']);

Y=CVF1INTRA;
Y_text='CV_F (nondim.)';
X=F1;
X_text='F (a.u.)';

h3=figure(3);


% plot data
for conditioni=1:4
    subplot(2,2,conditioni)
    plot(X(:,conditioni),Y(:,conditioni),'o','markersize',msize,'markerfacecolor',col(conditioni,:),'markeredgecolor',col(conditioni,:))
    hold on
    
    ylim([0.05 0.45]);
    xlim([0 6]);
    if conditioni==1
        set(gca,'xtick',[0 2 4],'xticklabel',[],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[0 "" 0.2 "" 0.4 "" 0.6]);
    elseif conditioni==2
        set(gca,'xtick',[0 2 4],'xticklabel',[],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[]);
    elseif conditioni==3
        set(gca,'xtick',[0 2 4 6],'xticklabel',[0 2 4 6],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[0 "" 0.2 "" 0.4 "" 0.6]);
    elseif conditioni==4
        set(gca,'xtick',[0 2 4 6],'xticklabel',[0 2 4 6],'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6],'yticklabel',[]);
	end
	set(gca,'fontsize',fsize);
    set(gca,'box','off') %Remove y right ticks    
    title([titles{:,conditioni}],'fontsize',fsize,'fontweight','normal')
end


% fit Linear Mixed Model (LMM)
Y=CVF1INTRA; X=1./F1;
YY=[Y(:,1);Y(:,2);Y(:,3);Y(:,4)];
XX=[X(:,1);X(:,2);X(:,3);X(:,4)];
CONf=[zeros(22,1);ones(22,1);zeros(22,1);ones(22,1)];
HIGH=[zeros(44,1);ones(44,1)];
subject=[[1:22]';[1:22]';[23:44]';[23:44]'];

tbl=table(subject,YY,CONf,HIGH,XX);
lme = fitlme(tbl,'YY ~ 1 + CONf + HIGH + XX + CONf:HIGH + XX:CONf +XX:HIGH+ (1|subject)');
lme_anova = anova(lme);

% slope and intercept for every condition
b(1)=lme.Coefficients.Estimate(1);
b(2)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2);
b(3)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(3);
b(4)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2)+lme.Coefficients.Estimate(3)+lme.Coefficients.Estimate(5);

m(1)=lme.Coefficients.Estimate(4);
m(2)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6);
m(3)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(7);
m(4)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6)+lme.Coefficients.Estimate(7);

% R_squared calculation
for i=1:4    
    y_mean=mean(Y(:,i));
    SS_tot=sum((Y(:,i)-y_mean).^2);
    SS_res=0;
    for j=1:length(Y(:,1))
        temp=(m(i)*X(j,i)+b(i)-Y(j,i))^2;
        SS_res=SS_res+temp;
    end
    RSquared_LMM=1-(SS_res/SS_tot);
    p=2; %number of parameters of the model (slope and intercept)
    L=length(Y(:,1)); %sample size
    RSquared_LMM_adj=1-(1-RSquared_LMM)*(L-1)/(L-p);
end


% confidence intervals of regression
% ref: https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals
n=size(X,1);
df=n-2; %degrees of freedom
t_score=tinv(probability,df);

% plot curve fit and calculate and plot confidence intervals
for conditioni=1:4
    %LMM fit
    subplot(2,2,conditioni);    
    slope=m(conditioni);
    intercept=b(conditioni);
    x=X(:,conditioni);
    y=Y(:,conditioni);
    xFit = linspace(min(x),max(x), 20);
    yFit = slope*xFit+intercept;
    yFit_onlyExperimentalValues=slope*x+intercept;
    hold on
    plot(1./xFit,yFit,'-','MarkerSize',1,'LineWidth',2,'Color','k')
    MSE=sum((y-yFit_onlyExperimentalValues).^2)/(n-2);
    SSx=sum((x-mean(x)).^2);
    DIF2=(xFit-mean(x)).^2;
    error=t_score*sqrt(MSE*(1/n+DIF2./SSx));
    upper=yFit+error;
    lower=yFit-error;
    hold on
    patch([1./xFit fliplr(1./xFit)], [lower fliplr(upper)], [0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none')
end



h3_ax = axes(h3,'visible','off');
h3_ax.Title.Visible = 'on';
h3_ax.XLabel.Visible = 'on';
h3_ax.YLabel.Visible = 'on';
xlabel(h3_ax,X_text);
ylabel(h3_ax,Y_text);
set(gca,'fontsize',fsize);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h3);



disp([newline '<strong>FIGURE 4: CVF vs F</strong>']);
disp(['R^2=' num2str(round(100*RSquared_LMM,4)/100)]);
disp(lme);
disp(lme_anova);
print(h3,'-dpdf','-bestfit','./figure5_CVF_vs_F.pdf');



%% Figure 6: MSJ vs F

% load data
load([data_dir 'F1.mat']);
load([data_dir 'MSJ.mat']);

% normalize and log
for conditioni=1:4
    temp_y=MSJ(:,conditioni);
    temp_x=F1(:,conditioni);
    Y(:,conditioni)=log(temp_y/max(temp_y));
    X(:,conditioni)=log(temp_x/max(temp_x));
end
Y_text='log (normalized MSJ) (nondim.)';
X_text='log (normalized F) (nondim.)';


h4=figure(4);


% plot data
for conditioni=1:4
    subplot(2,2,conditioni)
    plot(X(:,conditioni),Y(:,conditioni),'o','markersize',msize,'markerfacecolor',col(conditioni,:),'markeredgecolor',col(conditioni,:))

    %Spearman coeeficient and pValue
    x_noLog=F1(:,conditioni);
    y_noLog=MSJ(:,conditioni);
    [RHO(conditioni),PVAL] = corr(x_noLog,y_noLog,'Type','Spearman');
    hold on    
    xlim([-3.5 0.5])
    ylim([-10 0.7])
    
    if conditioni==1
        set(gca,'xtick',[-3 -2 -1 0],'xticklabel',[],'ytick',[-8 -4 0],'yticklabel',[-8 -4 0]);
    elseif conditioni==2
        set(gca,'xtick',[-3 -2 -1 0],'xticklabel',[],'ytick',[-8 -4 0],'yticklabel',[]);
    elseif conditioni==3
        set(gca,'xtick',[-3 -2 -1 0],'xticklabel',[-3 -2 -1 0],'ytick',[-8 -4 0],'yticklabel',[-8 -4 0]);
    elseif conditioni==4
        set(gca,'xtick',[-3 -2 -1 0],'xticklabel',[-3 -2 -1 0],'ytick',[-8 -4 0],'yticklabel',[]);
	end
	set(gca,'fontsize',fsize);
    set(gca,'box','off') %Remove y right ticks
    
    title([titles{:,conditioni}],'fontsize',fsize,'fontweight','normal')
               
end


% fit Linear Mixed Model (LMM)
YY=[Y(:,1);Y(:,2);Y(:,3);Y(:,4)];
XX=[X(:,1);X(:,2);X(:,3);X(:,4)];
CONf=[zeros(22,1);ones(22,1);zeros(22,1);ones(22,1)];
HIGH=[zeros(44,1);ones(44,1)];
subject=[[1:22]';[1:22]';[23:44]';[23:44]'];

tbl=table(subject,YY,CONf,HIGH,XX);
lme = fitlme(tbl,'YY ~ 1 + CONf + HIGH + XX + CONf:HIGH + XX:CONf +XX:HIGH+ (1|subject)');
lme_anova = anova(lme);

% slope and intercept for every condition
b(1)=lme.Coefficients.Estimate(1);
b(2)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2);
b(3)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(3);
b(4)=lme.Coefficients.Estimate(1)+lme.Coefficients.Estimate(2)+lme.Coefficients.Estimate(3)+lme.Coefficients.Estimate(5);

m(1)=lme.Coefficients.Estimate(4);
m(2)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6);
m(3)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(7);
m(4)=lme.Coefficients.Estimate(4)+lme.Coefficients.Estimate(6)+lme.Coefficients.Estimate(7);

% R_squared calculation
for i=1:4    
    y_mean=mean(Y(:,i));
    SS_tot=sum((Y(:,i)-y_mean).^2);
    SS_res=0;
    for j=1:length(Y(:,1))
        temp=(m(i)*X(j,i)+b(i)-Y(j,i))^2;
        SS_res=SS_res+temp;
    end
    RSquared_LMM=1-(SS_res/SS_tot);
    p=2; %number of parameters of the model (slope and intercept)
    L=length(Y(:,1)); %sample size
    RSquared_LMM_adj=1-(1-RSquared_LMM)*(L-1)/(L-p);
end


% confidence intervals of regression
% ref: https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals
n=size(X,1);
df=n-2; %degrees of freedom
t_score=tinv(probability,df);

% plot curve fit and calculate and plot confidence intervals
for conditioni=1:4
    %LMM fit
    subplot(2,2,conditioni)
    
    slope=m(conditioni);
    intercept=b(conditioni);
    x=X(:,conditioni);
    y=Y(:,conditioni);
    xFit = linspace(min(x),max(x), 20);
    yFit = slope*xFit+intercept;
    yFit_onlyExperimentalValues=slope*x+intercept;
    hold on
    plot(xFit,yFit,'-','MarkerSize',1,'LineWidth',2,'Color','k')
    MSE=sum((y-yFit_onlyExperimentalValues).^2)/(n-2);
    SSx=sum((x-mean(x)).^2);
    DIF2=(xFit-mean(x)).^2;
    error=t_score*sqrt(MSE*(1/n+DIF2./SSx));
    upper=yFit+error;
    lower=yFit-error;
    hold on
    patch([xFit fliplr(xFit)], [lower fliplr(upper)], [0.7 0.7 0.7],'facealpha',0.3,'edgecolor','none')
end



h4_ax = axes(h4,'visible','off');
h4_ax.Title.Visible = 'on';
h4_ax.XLabel.Visible = 'on';
h4_ax.YLabel.Visible = 'on';
xlabel(h4_ax,X_text);
ylabel(h4_ax,Y_text);
set(gca,'fontsize',fsize);

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h4);



disp([newline '<strong>FIGURE 5: MSJ vs F</strong>']);
disp(['R^2=' num2str(round(100*RSquared_LMM,4)/100)]);
disp(lme);
disp(lme_anova);
print(h4,'-dpdf','-bestfit','./figure6_MSJ_vs_F.pdf');



%% Supplementary Figure S1, force profiles


% load data
load([data_dir 'FORCE.mat']);
load([data_dir 'FMEAN.mat']);
load([data_dir 'MSJ.mat']);
load([data_dir 'data_timing_force.mat'],'NSujetos');


Y_text='Force (a.u.)';
X_text='Time (ms)';


h5=figure(5);

lwidth = 1;
fsize = 9;
% msize=3;
fig_size_cm = [20 14];


subplot(2,2,[1 2]);
hold on;
ptp = [];
conditioni = 3;
yeari = 1;
attentioni = 2;
subjecti = 4;
feedbacki = 1;
for perti = 1:2
    for triali = 1:6
        for tapi = 1:7
            tap = squeeze(force(yeari,attentioni,subjecti,feedbacki,perti,triali,tapi,:));
            ptp(1) = plot(1:size(tap,1),tap,'linewidth',1,'color',[0.8 0.8 0.8]);
% 			hold on;
        end
    end
end
set(gca,'ytick',[0 1 2 3 4],'yticklabel',[0 1 2 3 4],'xtick',[0 50 100 150],'fontsize',fsize)


subjecti = subjecti + (yeari-1)*NSujetos;
tap2 = squeeze(FMEAN(subjecti,:,conditioni));
ptp(2) = plot(1:size(tap2,2),tap2,'linewidth',2,'color',col(conditioni,:));

xlim([0 150])
ylim([0 4.5])
xlabel(X_text);
ylabel(Y_text);
legend(ptp,{'Individual taps','Mean tap'},'fontsize',fsize);
legend boxoff;
text(-10,4.3,'(A)','fontsize',fsize+2);





% MSJ
conditioni=1;

% Low MSJ
subplot(2,2,3);
subjecti=6;
plot(1:size(FMEAN,2),FMEAN(subjecti,:,conditioni),'k-','linewidth',2);
text(110,1.8,['(MSJ=' num2str(MSJ(subjecti,conditioni),'%.1E') ')'],'fontsize',fsize);
ylim([0 2]);
xlim([0 180]);
set(gca,'ytick',[1 2],'yticklabel',[1 2],'xtick',[0 50 100 150],'xticklabel',[0 50 100 150],'fontsize',fsize);
xlabel(X_text);
ylabel(Y_text);
text(-30,1.95,'(B)','fontsize',fsize+2);

% High MSJ
subplot(2,2,4);
subjecti=7;
plot(1:size(FMEAN,2),FMEAN(subjecti,:,conditioni),'k-','linewidth',2);
text(110,1.8,['(MSJ=' num2str(MSJ(subjecti,conditioni),'%.1E') ')'],'fontsize',fsize);
ylim([0 2]);
xlim([0 180]);
set(gca,'ytick',[1 2],'yticklabel',[1 2],'xtick',[0 50 100 150],'xticklabel',[0 50 100 150],'fontsize',fsize);
xlabel(X_text);
ylabel(Y_text);



set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h5);


print(h5,'-dpdf','-bestfit','./figureS1_forceprofiles.pdf');



%% POWER ANALYSIS, Supplementary Figure S2

% lag1 autocorrelation of asynchronies
% (see Semjen, Schulze & Vorberg 2000, figure 4G; Repp 2011, figure 8)
autocorr_lag1 = 0.4;

% data from Sasaki et al 2011
% SD of intertap interval
aux1_df = readtable([pow_dir 'Sasaki 2011 figure 3C.dat']); % ISI=250ms
aux2_df = readtable([pow_dir 'Sasaki 2011 figure 3D.dat']); % ISI=500ms
data_iti_Sasaki_df = [aux1_df; aux2_df];
% mean force
aux1_df = readtable([pow_dir 'Sasaki 2011 figure 5A.dat']); % ISI=250ms
aux2_df = readtable([pow_dir 'Sasaki 2011 figure 5B.dat']); % ISI=500ms
data_force_mean_Sasaki_df = [aux1_df; aux2_df];
% SD force
aux1_df = readtable([pow_dir 'Sasaki 2011 figure 5E.dat']); % ISI=250ms
aux2_df = readtable([pow_dir 'Sasaki 2011 figure 5F.dat']); % ISI=500ms
data_force_sd_Sasaki_df = [aux1_df; aux2_df];
% force CV
data_force_Sasaki_df = join(data_force_mean_Sasaki_df,data_force_sd_Sasaki_df);
data_force_Sasaki_df(:,'force_cv') = table(data_force_Sasaki_df(:,'force_sd').Variables./data_force_Sasaki_df(:,'force_mean').Variables);
% merge ITI and force
data_Sasaki_df = join(data_iti_Sasaki_df,data_force_Sasaki_df);
data_Sasaki_df(:,'asyn_sd') = table(data_Sasaki_df(:,'iti_sd').Variables/sqrt(2*(1-autocorr_lag1)));

% data from Inui et al 2002
% cv of force
data_force_Inui_df = readtable([pow_dir 'Inui 2002 figure 5A.dat']);
% ITI
aux1_df = readtable([pow_dir 'Inui 2002 figure 2B.dat']); % mean intertap interval
aux2_df = readtable([pow_dir 'Inui 2002 figure 5B.dat']); % cv of intertap interval
data_iti_Inui_df = join(aux1_df,aux2_df);
data_iti_Inui_df(:,'iti_sd') = table(data_iti_Inui_df(:,'iti_mean').Variables.*data_iti_Inui_df(:,'iti_cv').Variables);
data_iti_Inui_df(:,'asyn_sd') = table(data_iti_Inui_df(:,'iti_sd').Variables/sqrt(2*(1-autocorr_lag1)));
% merge ITI and force
data_Inui_df = join(data_iti_Inui_df,data_force_Inui_df);
% remove tap=4 due to different experimental condition
data_Inui_df(4,:) = [];

% merge all data
data_df = [data_Sasaki_df(:,{'force_cv','asyn_sd'}); data_Inui_df(:,{'force_cv','asyn_sd'})];

% estimate correlation
force_cv = data_df(:,'force_cv').Variables;
asyn_sd = data_df(:,'asyn_sd').Variables;
corr_SDA_CVF = corr(asyn_sd,force_cv, 'type','pearson');
% linear regression
model = fitlm(data_df(:,{'force_cv','asyn_sd'}));

% estimate number of subjects
% (Suresh & Chandrashekara, 2012; section “Sample size estimation with correlation coefficient”)
z_onetailed = 1.64;
z_pow80 = 0.84;
N_estim = ceil(4*(z_onetailed+z_pow80)^2/log((1+corr_SDA_CVF)/(1-corr_SDA_CVF)) + 3);

disp(['Correlation coefficient SD_A vs CV_F from literature= ' num2str(corr_SDA_CVF)]);
disp(['Estimated number of participants: ' num2str(N_estim)]);



% plot data and regression
h6=figure(6);

fsize=8;
msize = 10;
fig_size_cm=[12 8];

ptp = plot(model);
ptp(1).Marker = '.';
ptp(1).MarkerSize = msize;
ptp(1).Color = col(3,:);
ptp(1).MarkerFaceColor = 'm';
ptp(2).Color = 'k';
ptp(2).LineWidth = 1.25;
ptp(3).Color = 'k';
ptp(3).LineWidth = 1;
ptp(3).LineStyle = ':';
ptp(4).Color = 'k';
ptp(4).LineWidth = 1;
ptp(4).LineStyle = ':';

set(gca,'fontsize', fsize);
xlim([0.13 0.32]);
ylim([0 70]);
xlabel('CV_F (nondim.)','fontsize',fsize,'interpreter','tex');
ylabel('SD_A (ms)','fontsize',fsize,'interpreter','tex');
title('');
legend(ptp(1),['Data from' newline 'literature'],'location','northwest');
legend boxoff;


set(gcf,'units','centimeters','position',[0,0,fig_size_cm]);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

tightfig(h6);


print(h6,'-dpdf','-bestfit','./figureS2_previousdata.pdf');




%%

