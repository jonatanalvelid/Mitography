%%% MULTIPLE HISTOGRAM PLOTTING
% Dataset: \\X:\TestaLab\Mitography\MitoSOX-MitographyAnalysis\mitoData-RL-200330.mat

colors = lines(2);
gray = [0.6 0.6 0.6];
lightGray = [0.7 0.7 0.7];
darkGray = [0.2 0.2 0.2];

% Take all mitos
% h1var = mitoWidth(mitoMitoSOXparam==1);  % All MitoSOX+ mito
% h2var = mitoWidth(mitoMitoSOXparam==0);  % All MitoSOX- mito
% h3var = mitoArea(mitoMitoSOXparam==1);  % All MitoSOX+ mito
% h4var = mitoArea(mitoMitoSOXparam==0);  % All MitoSOX- mito
% h5var = mitoLength(mitoMitoSOXparam==1);  % All MitoSOX+ mito
% h6var = mitoLength(mitoMitoSOXparam==0);  % All MitoSOX- mito

% Only take small mitos
areathresh = 0.086;
mitoWidthsmall = mitoWidth(mitoArea<areathresh);
mitoWidthbig = mitoWidth(mitoArea>areathresh);
mitoLengthsmall = mitoLength(mitoArea<areathresh);
mitoLengthbig = mitoLength(mitoArea>areathresh);
mitoAreasmall = mitoArea(mitoArea<areathresh);
mitoAreabig = mitoArea(mitoArea>areathresh);
mitoMitoSOXparamsmall = mitoMitoSOXparam(mitoArea<areathresh);
mitoMitoSOXparambig = mitoMitoSOXparam(mitoArea>areathresh);

h1var = mitoWidthsmall(mitoMitoSOXparamsmall==1);  % Small MitoSOX+ mito
h2var = mitoWidthsmall(mitoMitoSOXparamsmall==0);  % Small MitoSOX- mito
h3var = mitoAreasmall(mitoMitoSOXparamsmall==1);  % Small MitoSOX+ mito
h4var = mitoAreasmall(mitoMitoSOXparamsmall==0);  % Small MitoSOX- mito
h5var = mitoLengthsmall(mitoMitoSOXparamsmall==1);  % Small MitoSOX+ mito
h6var = mitoLengthsmall(mitoMitoSOXparamsmall==0);  % Small MitoSOX- mito
% h3var = h1var.*h5var;  % Small MitoSOX+ mito
% h4var = h2var.*h6var;  % Small MitoSOX- mito

boundlow = [0, 0, 0, 0];
stepwidth = [0.015, 0.005, 0.04];
boundup = [0.4, 0.1, 20];
xlimlow = boundlow;
xlimup = boundup;
ylimup = [0.4, 0.3, 0.3];

fontsize = 12;
opacity = 0.5;

% legendtext1 = 'All mito';
legendtext1 = sprintf('MitoSOX+');
legendtext2 = sprintf('MitoSOX-');

xlabeltext1 = 'Mitochondria width [um]';
xlabeltext2 = 'Mitochondria area [um^2]';
xlabeltext3 = 'Mitochondria length [um]';


%%%,'FaceAlpha',opacity %%% If you want different opacity
mitowidfig = figure('rend','painters','pos',[100 100 300 300]);
n = 1;
h1 = histogram(h1var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',lightGray);
hold on
h2 = histogram(h2var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',darkGray);
xlim([xlimlow(n) xlimup(n)])
xlabel(xlabeltext1)
ylim([0 ylimup(n)])
ylabel('Norm. frequency')
%title(strcat(titletext1,', N=',num2str(length(h1var))));
legend(legendtext1,legendtext2);
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})
yticks([0:ylimup(n)/12:ylimup(n)])
yticklabels({0,'','',ylimup(n)/4,'','',ylimup(n)/2,'','',3*ylimup(n)/4,'','',ylimup(n)})
disp(mean(h1var))
disp(mean(h2var))
disp(median(h1var))
disp(median(h2var))
disp(' ')

mitoareafig = figure('rend','painters','pos',[500 100 300 300]);
n = 2;
h3 = histogram(h3var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',lightGray);
hold on
h4 = histogram(h4var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',darkGray);
xlim([xlimlow(n) xlimup(n)])
xlabel(xlabeltext2)
ylim([0 ylimup(n)])
ylabel('Norm. frequency')
%title(strcat(titletext2,', N=',num2str(length(h3var))));
legend(legendtext1,legendtext2);
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})
yticks([0:ylimup(n)/12:ylimup(n)])
yticklabels({0,'','',ylimup(n)/4,'','',ylimup(n)/2,'','',3*ylimup(n)/4,'','',ylimup(n)})
disp(mean(h3var))
disp(mean(h4var))
disp(median(h3var))
disp(median(h4var))
disp(' ')

mitolenfig = figure('rend','painters','pos',[900 100 300 300]);
n = 3;
h5 = histogram(h5var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',lightGray);
hold on
h6 = histogram(h6var,boundlow(n):stepwidth(n):boundup(n),'Normalization','probability','FaceColor',darkGray);
xlim([xlimlow(n) xlimup(n)])
xlabel(xlabeltext3)
ylim([0 ylimup(n)])
ylabel('Norm. frequency')
%title(strcat(titletext3,', N=',num2str(length(h5var))));
legend(legendtext1,legendtext2);
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})
yticks([0:ylimup(n)/12:ylimup(n)])
yticklabels({0,'','',ylimup(n)/4,'','',ylimup(n)/2,'','',3*ylimup(n)/4,'','',ylimup(n)})
disp(mean(h5var))
disp(mean(h6var))
disp(median(h5var))
disp(median(h6var))
disp(' ')

%{
% mitoarealogfig = figure('rend','painters','pos',[1700 100 300 300]);
% n = 5; nbins = 30;
% [~,binedges] = histcounts(log10(h9var),nbins);
% histogram(h9var,10.^binedges,'Normalization','probability')
% hold on
% [~,binedges] = histcounts(log10(h10var),nbins);
% histogram(h10var,10.^binedges,'Normalization','probability')
% % [~,binedges] = histcounts(log10(h11var),nbins);
% % histogram(h11var,10.^binedges,'Normalization','probability')
% % [~,binedges] = histcounts(log10(h12var),nbins);
% % histogram(h12var,10.^binedges,'Normalization','probability')
% set(gca, 'xscale','log')
% ylabel('Norm. frequency')
% xlabel(xlabeltext2)
% %title(strcat(titletext3,', N=',num2str(length(h5var))));
% % legend(legendtext1,legendtext2);
% legend('Raw','RL');
% set(gca,'FontSize',fontsize)
% set(gca,'TickDir','out');


% Box plots
hdouble = [h1var;h2var];
groupings = [zeros(1,length(h1var)),ones(1,length(h2var))];
widthboxplots = figure('rend','painters','pos',[100 100 300 150]);
n = 1;
h9 = boxplot(hdouble,groupings,'orientation','horizontal','BoxStyle','outline','Colors',colors,'OutlierSize',4,'Symbol','k+','Widths',0.8);
xlim([xlimlow(n) xlimup(n)])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'yticklabel',[])
set(gca,'ytick',[])
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})

hdouble = [h3var;h4var];
groupings = [zeros(1,length(h3var)),ones(1,length(h4var))];
areaboxplots = figure('rend','painters','pos',[100 100 300 150]);
n = 2;
h10 = boxplot(hdouble,groupings,'orientation','horizontal','BoxStyle','outline','Colors',colors,'OutlierSize',4,'Symbol','k+','Widths',0.8);
xlim([xlimlow(n) xlimup(n)])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'yticklabel',[])
set(gca,'ytick',[])
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})

hdouble = [h5var;h6var];
groupings = [zeros(1,length(h5var)),ones(1,length(h6var))];
lengthboxplots = figure('rend','painters','pos',[100 100 300 150]);
n = 3;
h11 = boxplot(hdouble,groupings,'orientation','horizontal','BoxStyle','outline','Colors',colors,'OutlierSize',4,'Symbol','k+','Widths',0.8);
xlim([xlimlow(n) xlimup(n)])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'yticklabel',[])
set(gca,'ytick',[])
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})

hdouble = [h7var;h8var];
groupings = [zeros(1,length(h7var)),ones(1,length(h8var))];
arboxplots = figure('rend','painters','pos',[100 100 300 150]);
n = 4;
h12 = boxplot(hdouble,groupings,'orientation','horizontal','BoxStyle','outline','Colors',colors,'OutlierSize',4,'Symbol','k+','Widths',0.8);
xlim([xlimlow(n) xlimup(n)])
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
set(gca,'yticklabel',[])
set(gca,'ytick',[])
xticks([xlimlow(n):(xlimup(n)-xlimlow(n))/12:xlimup(n)])
xticklabels({xlimlow(n),'','',(xlimup(n)-xlimlow(n))/4,'','',(xlimup(n)-xlimlow(n))/2,'','',3*(xlimup(n)-xlimlow(n))/4,'','',xlimup(n)})
%}