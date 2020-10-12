%%%
% Mitography - Turnover analysis - data plotting
% Plotting the data from the turnover analysis.
%----------------------------
% Version: 201012
% Last updated features: New script
%
% @jonatanalvelid
%%%

clear

% Add function folder to filepath, so that those functions can be read.
functionFolder = fileparts(which('findFunctionFolders.m'));
addpath(genpath(functionFolder));

masterFolderPath = strcat(uigetdir('T:\Mitography'),'\');
fileList = dir(fullfile(masterFolderPath, '*turnoveranalysis.txt'));
filenumbers = [];
for i = 1:length(fileList)
    filenumbers(i) = str2num(fileList(i).name(1:3));
end
lastFileNumber = max(filenumbers);

mitodatatemp = [];
mitodata = [];
somadist = [];
signratio = [];
area = [];

filenameallMito = '_turnoveranalysis.txt';
fileNumbers = 1:lastFileNumber;

%%% GATHER DATA

for fileNum = fileNumbers
    filepathmito = strFilepath2(fileNum,filenameallMito,masterFolderPath);
    
    try
        % Read the mito and line profile data
        datamito = dlmread(filepathmito,'',0,0);
        [num,params] = size(datamito);
        
        % Get somadist and number of nucleotides for all mito in image
        somabooltemp = datamito(:,params-3);
        somadisttemp = datamito(:,params-2);
        tmrsignaltemp = datamito(:,params-1);
        sirsignaltemp = datamito(:,params);
        signratiotemp = sirsignaltemp./tmrsignaltemp;
        areatemp = datamito(:,4);
        
        % Add somadist = 0 if somabool = 1
        somadisttemp(somabooltemp==1) = 0;
        
        % Remove all other NaN somadist mitos
        idxrem = isnan(somadisttemp);
        somadisttemp(idxrem) = [];
        signratiotemp(idxrem) = [];
        areatemp(idxrem) = [];
        
        % Add to list of all mitochondria in all images
        somadist = vertcat(somadist,somadisttemp);
        signratio = vertcat(signratio,signratiotemp);
        area = vertcat(area,areatemp);
        
    catch err
        fprintf(1,'Line of error:\n%i\n',err.stack(end).line);
        fprintf(1,'The identifier was:\n%s\n',err.identifier);
        fprintf(1,'There was an error! The message was:\n%s\n',err.message);
    end 
end

areathlo = 0.03;
areathhi = 0.15;
%areathlo = 0.03;
%areathhi = 5;

mitodatatemp(:,1) = somadist(area<areathhi & area>areathlo);
mitodatatemp(:,2) = signratio(area<areathhi & area>areathlo);
mitodatatemp(:,3) = area(area<areathhi & area>areathlo);

mitodata = vertcat(mitodata,mitodatatemp);
xdata = mitodata(:,1);
ydata = mitodata(:,2);
ydata(ydata>10) = 10;
ydata(ydata<0.01) = 0.01;

%%% PLOTTING CODE

colors = lines(2);
gray = [0.6 0.6 0.6];
lightGray = [0.7 0.7 0.7];
darkGray = [0.3 0.3 0.3];

xlimlow1 = 0;
xlimup1 = 130;
ylimlow1 = 0.01;
ylimup1 = 10;

fontsize = 14;
opacity = 0.5;

ylabeltext1 = 'Signal ratio, SiR/TMR (a.u.)';
xlabeltext1 = 'Distance to soma (um)';

TFAMfig = figure('rend','painters','pos',[100 100 800 400]);
h1 = scatter(xdata,ydata);
xlim([xlimlow1 xlimup1])
ylim([ylimlow1 ylimup1])
xlabel(xlabeltext1)
ylabel(ylabeltext1)
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
xticks([xlimlow1:(xlimup1-xlimlow1)/12:xlimup1])
xticklabels({xlimlow1,'','',(xlimup1-xlimlow1)/4,'','',(xlimup1-xlimlow1)/2,'','',3*(xlimup1-xlimlow1)/4,'','',xlimup1})
set(gca, 'YScale', 'log')
yticks([0.01 0.1 1 10 100])
yticklabels({'0.01','green','yellow','red','100'})

[~,edges] = histcounts(log10(ydata));

figure(2)
histogram(ydata(xdata<10),10.^edges,'Normalization','probability')
hold on
histogram(ydata(xdata>10 & xdata<50),10.^edges,'Normalization','probability')
histogram(ydata(xdata>50),10.^edges,'Normalization','probability')
set(gca, 'xscale','log')
xlim([ylimlow1 ylimup1])