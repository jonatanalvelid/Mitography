clear all

dirread = strcat('E:\PhD\data_analysis\temp_pex','\');
datalist = dir(fullfile(dirread,'*.mat'));

mdvoxpn_area = [];
mdvoxpp_area = [];
mdvoxp_rat = [];
posratios = [];
% PEX exp 2
%expnum = [1 1 1 1];
%cellnum = [1 2 3 4];
% PEX - all exp
expnum = [1 1 1 1 1 2 2 3 3 3 3];
cellnum = [1 2 3 4 5 1 2 1 2 3 4];

area_pooled = struct;
for i=1:max(expnum)
    area_pooled(i).areasn = [];
    area_pooled(i).areasp = [];
end

figure()
for i=1:numel(datalist)
    disp(' ')
    disp(datalist(i).name)
    
    exp = expnum(i);
    cell = cellnum(i);
    
    load(strcat(datalist(i).folder,'\',datalist(i).name));
    
    idx_n = mitoinfo.area<0.086 & mitoinfo.pexparam == 0;
    idx_p = mitoinfo.area<0.086 & mitoinfo.pexparam == 1;
    meanarea_n = mean(mitoinfo.area(idx_n));
    meanarea_p = mean(mitoinfo.area(idx_p));
    mdvoxpn_area = [mdvoxpn_area meanarea_n];
    mdvoxpp_area = [mdvoxpp_area meanarea_p];
    mdvoxp_rat = [mdvoxp_rat sum(idx_p)/(sum(idx_p)+sum(idx_n))];
    fprintf('Mean area, -:  %f \n', meanarea_n)
    fprintf('Mean area, +:  %f \n', meanarea_p)
    
    n_n = length(mitoinfo.area(idx_n));
    n_p = length(mitoinfo.area(idx_p));
    pos_ratio = n_p/(n_p+n_n);
    posratios = [posratios pos_ratio];
    fprintf('+ ratio of MDVs:  %.3f \n', pos_ratio)
    
    area_pooled(exp).areasn = [area_pooled(exp).areasn; mitoinfo.area(idx_n)];
    area_pooled(exp).areasp = [area_pooled(exp).areasp; mitoinfo.area(idx_p)];
    try
        [~,p] = kstest2(mitoinfo.area(idx_n),mitoinfo.area(idx_p));
        fprintf('KS-test (+/-): %f \n', p)
    catch
    end
    
    subplot(max(expnum),max(cellnum),(exp-1)*max(cellnum)+cell);
    histogram(mitoinfo.area(idx_n),[0:0.01:0.1],'FaceColor',[0.3 0.3 0.3],'Normalization','probability')
    hold on
    histogram(mitoinfo.area(idx_p),[0:0.01:0.1],'FaceColor',[0 0.7 0],'Normalization','probability')
    xlim([0 0.1])
    xlabel('Area (um^2)')
    %ylim([0 1])
    ylabel('Probability (arb.u.)')
    if exp == 1
        title(sprintf('Exp: %d - AA, Cell: %d,\n N+=%d, N-=%d, p=%.3f',exp,cell,length(mitoinfo.area(idx_p)),length(mitoinfo.area(idx_n)),p))
    else
        title(sprintf('Exp: %d - CTR, Cell: %d,\n N+=%d, N-=%d, p=%.3f',exp-1,cell,length(mitoinfo.area(idx_p)),length(mitoinfo.area(idx_n)),p))
    end
end

figure()
for i=1:numel(datalist)
    disp(' ')
    disp(datalist(i).name)
    
    exp = expnum(i);
    cell = cellnum(i);
    
    load(strcat(datalist(i).folder,'\',datalist(i).name));
    
    idx_n = mitoinfo.area<0.086 & mitoinfo.pexparam == 0;
    idx_p = mitoinfo.area<0.086 & mitoinfo.pexparam == 1;
    meanarea_n = mean(mitoinfo.area(idx_n));
    meanarea_p = mean(mitoinfo.area(idx_p));
    mdvoxpn_area = [mdvoxpn_area meanarea_n];
    mdvoxpp_area = [mdvoxpp_area meanarea_p];
    mdvoxp_rat = [mdvoxp_rat sum(idx_p)/(sum(idx_p)+sum(idx_n))];
    fprintf('Mean area, -:  %f \n', meanarea_n)
    fprintf('Mean area, +:  %f \n', meanarea_p)
    
    n_n = length(mitoinfo.area(idx_n));
    n_p = length(mitoinfo.area(idx_p));
    try
        if n_n ~= 0 && n_p ~= 0
            subplot(max(expnum),max(cellnum),(exp-1)*max(cellnum)+cell);
            hold on
            labels = {'pex+', 'pex-'};
            groups = [ones(n_p,1);2*ones(n_n,1)];
            areap = mitoinfo.area(idx_p);
            arean = mitoinfo.area(idx_n);
            boxplot([areap;arean],groups,'Labels',labels)
            s1 = scatter(ones(size(areap)),areap,10,'k.','jitter','on','jitterAmount',0.1);
            s2 = scatter(2*ones(size(arean)),arean,10,'k.','jitter','on','jitterAmount',0.1);
            [~,p] = kstest2(areap,arean);
        end
    catch
    end
    ylabel('Area (um^2)')
    ylim([0 0.1])
    if exp == 1
        title(sprintf('Exp: %d - AA, Cell: %d,\n N+=%d, N-=%d, p=%.3f',exp,cell,length(areap),length(arean),p))
    else
        title(sprintf('Exp: %d - CTR, Cell: %d,\n N+=%d, N-=%d, p=%.3f',exp-1,cell,length(areap),length(arean),p))
    end
end

figure()
subplotpos = [1 3 4];
exps = [1 2 3];
titles = {'Exp: %d - AA,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f'};
for i=1:max(expnum)
    subplot(2,2,subplotpos(i));
    hold on
    histogram(area_pooled(i).areasn,[0:0.005:0.1],'FaceColor',[0.3 0.3 0.3],'Normalization','probability')
    histogram(area_pooled(i).areasp,[0:0.005:0.1],'FaceColor',[0 0.7 0],'Normalization','probability')
    
    [~,p] = kstest2(area_pooled(i).areasn,area_pooled(i).areasp);
    fprintf('KS-test (+/-): %f \n', p)
    
    xlim([0 0.1])
    xlabel('Area (um^2)')
    ylim([0 0.4])
    ylabel('Probability (arb.u.)')
    title(sprintf(titles{i},exps(i),length(area_pooled(i).areasp),length(area_pooled(i).areasn),p))
end

figure()
subplotpos = [1 3 4];
exps = [1 2 3];
titles = {'Exp: %d - AA,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f','Exp: %d - CTR,\n N+=%d, N-=%d, p=%.3f'};
for i=1:max(expnum)
    subplot(2,2,subplotpos(i));
    hold on
    labels = {'pex+', 'pex-'};
    areap = area_pooled(i).areasp;
    arean = area_pooled(i).areasn;
    n_n = length(arean);
    n_p = length(areap);
    groups = [ones(n_p,1);2*ones(n_n,1)];
    boxplot([areap;arean],groups,'Labels',labels)
    s1 = scatter(ones(size(areap)),areap,30,'k.','jitter','on','jitterAmount',0.1);
    s2 = scatter(2*ones(size(arean)),arean,30,'k.','jitter','on','jitterAmount',0.1);
            
    [~,p] = kstest2(areap,arean);
    fprintf('KS-test (+/-): %f \n', p)
    ylabel('Area (um^2)')
    ylim([0 0.1])
    title(sprintf(titles{i},exps(i),length(areap),length(arean),p))
end

disp(' ')
for n=1:max(expnum)
    for m=1:max(expnum)
        if n ~= m && m>n
            [~,p] = kstest2(area_pooled(n).areasp,area_pooled(m).areasp);
            fprintf('KS-test, %d v %d, (+/+): %f \n', n, m, p)
            [~,p] = kstest2(area_pooled(n).areasn,area_pooled(m).areasn);
            fprintf('KS-test, %d v %d, (-/-): %f \n', n, m, p)
        end
    end
end

% Box plots with jittered scatter of PEX+ ratios
figure()
hold on
labels = {'Exp1 - AA', 'Exp1 - CTR', 'Exp2 - CTR'};
s1 = scatter(ones(size(posratios(expnum==1))),posratios(expnum==1),40,'kx','jitter','on','jitterAmount',0.09);
s2 = scatter(2*ones(size(posratios(expnum==2))),posratios(expnum==2),40,'kx','jitter','on','jitterAmount',0.09);
s3 = scatter(3*ones(size(posratios(expnum==3))),posratios(expnum==3),40,'kx','jitter','on','jitterAmount',0.09);
boxplot(posratios,expnum,'Labels',labels)
ylabel('Ratio of PEX+ MDVs')
ylim([0 1])