%%%
% Mitography - TFAM analysis - summarize data
% Summarize the data and group based on axon distance.
%----------------------------
% Version: 191217
% Last updated features: New script
%
% @jonatanalvelid
%%%

groupsize = 20;

AxonDistance = mitodata(:,1);
NumNucleotides = mitodata(:,2);
AxonDistGroup = ceil(AxonDistance/groupsize)*groupsize;
AxonDistGroup(AxonDistGroup>300) = 300;  % The group for all distant mitos

mitodatatable = table(AxonDistance,AxonDistGroup,NumNucleotides);

stats = grpstats(mitodatatable,'AxonDistGroup',{'mean','sem','std'},'DataVars','NumNucleotides');

meannucl = table2array(stats(:,'mean_NumNucleotides'));
stdnucl = table2array(stats(:,'std_NumNucleotides'));
axondist = table2array(stats(:,'AxonDistGroup'));

% Check percentage of group elements that are 0
groups = groupsize:groupsize:300;
for i = 1:length(groups)
    numempty(i) = sum(NumNucleotides(AxonDistGroup==groups(i))==0);
    numtot(i) = length(NumNucleotides(AxonDistGroup==groups(i)));
    disp(groups(i))
    disp(numempty(i)/numtot(i))
end

scatter(groups, numempty./numtot)
ylim([0 1])
xlabel('Distance along axon [µm]')
ylabel('Ratio of mito with 0 nucleotides')
set(gca,'FontSize',14)


%%% PLOTTING CODE

colors = lines(2);
gray = [0.6 0.6 0.6];
lightGray = [0.7 0.7 0.7];
darkGray = [0.3 0.3 0.3];

xlimlow1 = 0;
xlimup1 = 300;
ylimup1 = 2;

fontsize = 14;

ylabeltext1 = 'Mean # of nucleotides';
xlabeltext1 = 'Distance along axon [µm]';

TFAMfig = figure('rend','painters','pos',[100 100 600 400]);
h1 = errorbar(axondist,meannucl,stdnucl);
xlim([xlimlow1 xlimup1])
ylim([0 ylimup1])
xlabel(xlabeltext1)
ylabel(ylabeltext1)
set(gca,'FontSize',fontsize)
set(gca,'TickDir','out');
% xticks([xlimlow1:(xlimup1-xlimlow1)/12:xlimup1])
% xticklabels({xlimlow1,'','',(xlimup1-xlimlow1)/4,'','',(xlimup1-xlimlow1)/2,'','',3*(xlimup1-xlimlow1)/4,'','',xlimup1})
% yticks([0:ylimup1/12:ylimup1])
% yticklabels({0,'','',ylimup1/4,'','',ylimup1/2,'','',3*ylimup1/4,'','',ylimup1})
