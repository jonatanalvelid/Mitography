function [filamentWidthSNZ, mitoDisplacementEdgeS, mitoActinDistS, mitoActinDistD2S, filamentWidthNSNZ, mitoDisplacementEdgeNS,...
    mitoActinDistNS, mitoActinDistD2NS, mitoWidthS, mitoAreaS, mitoLengthS, nonzeroS, mitoWidthNS, mitoAreaNS, mitoLengthNS, nonzeroNS, mitoNonStripesDistS,...
    mitoNonStripesDistNS, mitoThreePeaksNonStripesDist, mitoWidthRatioUMS, mitoWidthRatioBMS, mitoWidthRatioUBS, mitoWidthRatioUMNS,...
    mitoWidthRatioBMNS, mitoWidthRatioUBNS, mitoPatchOverlapS, mitoPatchOverlapNS, filamentWidthPNZ, mitoDisplacementEdgeP, mitoActinDistP,...
    mitoActinDistD2P, filamentWidthNPNZ, mitoDisplacementEdgeNP, mitoActinDistNP, mitoActinDistD2NP, mitoWidthP, mitoAreaP, mitoLengthP,...
    nonzeroP, mitoWidthNP, mitoAreaNP, mitoLengthNP, nonzeroNP, mitoPatchOverlapP, mitoPatchOverlapNP, mitoNonStripesDistP, mitoNonStripesDistNP, paramName] = updateData180612(doubleTol, doubleTol2, singleTol)

    masterFolderPath = strcat(uigetdir('D:\Data analysis\Mitography'),'\');
    filenameparam = 'ImageJAnalysisParameters.txt';
    filepathparam = strcat(masterFolderPath,filenameparam);   
    try
        dataparam = dlmread(filepathparam,'',1,1);
        mitoLen = dataparam(1,3);
        actinLen = dataparam(1,1);
    catch err
        mitoLen = 0;
        actinLen = 0;
    end
    onlyPlottingGUI(actinLen, mitoLen, 31, 0, 1000, 0.2, 2, 'Extra parameter')
    uiwait
    load('checkboxes.mat');
    load('edittextvals.mat');

    allCheckboxes = cell2mat(checks);
    paramPlots = allCheckboxes(2);
    doubleVals = str2double(vals(1:2));
    mitoLineProfLength = doubleVals(1);
    actinLineProfLength = doubleVals(2);
    intVals = str2double(vals(3:7));
    lastFileNumber = round(intVals(1));
    paramCol = round(intVals(2));
    mitosPerFile = round(intVals(3));
    mitoWidthRatioThresholdLow = intVals(4);
    mitoWidthRatioThresholdHigh = intVals(5);
    paramName = char(vals(8));

    %mitosPerFile = 1000;

    fileNumbers = 1:lastFileNumber;
    %c = 0:0.001:1;

    filenameAnalysis = '_MitoAnalysisFull.txt';
    filenameThreePeaks = '_MitoThreePeaks.txt';

    %%% DO ALL THE PLOTTING BELOW

    mitoThreePeaksWidthOneFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoThreePeaksWidthTwoFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoThreePeaksCenterOneFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoThreePeaksCenterTwoFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoThreePeaksNonStripesDistFiles = zeros(mitosPerFile,2,lastFileNumber);

    %%% TAKE CARE OF THE THREE PEAKS MITOCHONDRIA DATA

    for fileNum = fileNumbers
        if fileNum < 10
            filename = strcat('Image_00',int2str(fileNum),filenameThreePeaks);
        else
            filename = strcat('Image_0',int2str(fileNum),filenameThreePeaks);
        end
        filepath = strcat(masterFolderPath,filename);

        try
            data = dlmread(filepath);
            widthMitoOne = data(1:end,2);
            widthMitoTwo = data(1:end,4);
            centerMitoOne = data(1:end,3);
            centerMitoTwo = data(1:end,5);
            nonstripesdistMito = data(1:end,7);
            for i=1:length(widthMitoOne)
                mitoThreePeaksWidthOneFiles(i,1,fileNum) = i;
                mitoThreePeaksWidthOneFiles(i,2,fileNum) = widthMitoOne(i);
                mitoThreePeaksWidthTwoFiles(i,1,fileNum) = i;
                mitoThreePeaksWidthTwoFiles(i,2,fileNum) = widthMitoTwo(i);
                mitoThreePeaksCenterOneFiles(i,1,fileNum) = i;
                mitoThreePeaksCenterOneFiles(i,2,fileNum) = centerMitoOne(i);
                mitoThreePeaksCenterTwoFiles(i,1,fileNum) = i;
                mitoThreePeaksCenterTwoFiles(i,2,fileNum) = centerMitoTwo(i);
                mitoThreePeaksNonStripesDistFiles(i,1,fileNum) = i;
                mitoThreePeaksNonStripesDistFiles(i,2,fileNum) = nonstripesdistMito(i);
            end
        catch err
        end  
    end

    %%% TAKE CARE OF ALL MITOCHONDRIA DATA, DIVIDE FOR AIS AND STRIPES, AND
    %%% THREE PEAKS

    mitoWidthFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoLengthFiles = zeros(mitosPerFile,2,lastFileNumber);
    %mitoLengthSkelFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoAreaFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoCenterFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoNonStripesDistFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoStripesParam = zeros(mitosPerFile,2,lastFileNumber);
    mitoAISParam = zeros(mitosPerFile,2,lastFileNumber);
    mitoSmallParam = zeros(mitosPerFile,2,lastFileNumber);
    mitoBiggerParam = zeros(mitosPerFile,2,lastFileNumber);
    mitoRatioParamUMFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoRatioParamBMFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoRatioParamUBFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentWidthFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentCenterFiles = zeros(mitosPerFile,2,lastFileNumber);
    extraParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentDoublePeakParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentSinglePeakParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoSmallAreaParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoActinDistFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoActinDistNegParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoActinDistD2Files = zeros(mitosPerFile,2,lastFileNumber);
    mitoActinDistD2CloseParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoPatchOverlapFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentWidthMidFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentCenterMidFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentWidthUpFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentCenterUpFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentWidthBoFiles = zeros(mitosPerFile,2,lastFileNumber);
    filamentCenterBoFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoWidthMidFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoCenterMidFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoWidthUpFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoCenterUpFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoWidthBoFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoCenterBoFiles = zeros(mitosPerFile,2,lastFileNumber);


    for fileNum = fileNumbers
        if fileNum < 10
            filename = strcat('Image_00',int2str(fileNum),filenameAnalysis);
        else
            filename = strcat('Image_0',int2str(fileNum),filenameAnalysis);
        end
        filepath = strcat(masterFolderPath,filename);

        try
            data = dlmread(filepath);
            areaMito = data(1:end,4);
            lengthMito = data(1:end,5); %Ellipsoidal fit mitochondria length (major axis)
            lengthSkelMito = data(1:end,7); %Skeleton mitochondria length (skeleton closest to the mitochondria centroid)
            widthMito = data(1:end,8);
            centerMito = data(1:end,9);
            nonstripesdistMito = data(1:end,10); %Really patch-distance now, instead of non-stripes dist.
            actinwidth = data(1:end,11);
            actincenter = data(1:end,12);
            stripesparam = data(1:end,13);
            aisparam = data(1:end,14);
            smallparam = data(1:end,15);
            mitobiggerparam = data(1:end,16);
            mitoratioparamum = data(1:end,17);
            mitoratioparambm = data(1:end,18);
            mitoratioparamub = data(1:end,19);
            smallareaparam = data(1:end,20);
            actindoublepeakparam = data(1:end,21);
            actinsinglepeakparam = data(1:end,22);
            mitoactindist = data(1:end,24);
            mitoactindistnegparam = data(1:end,25);
            mitopatchoverlap = data(1:end,27);
            mitopatchdist = data(1:end,39);
            mitoactindistD2 = data(1:end,103);
            mitoactindistD2closeparam = data(1:end,104);
            
            widthMidMito = data(1:end,45);
            widthMidMito = zeros(size(widthMidMito));
            widthUpperMito = data(1:end,31);
            widthUpperMito = zeros(size(widthUpperMito));
            widthBottomMito = data(1:end,32);
            widthBottomMito = zeros(size(widthBottomMito));
            centerMidMito = data(1:end,46);
            centerMidMito = zeros(size(centerMidMito));
            centerUpperMito = data(1:end,61);
            centerUpperMito = zeros(size(centerUpperMito));
            centerBottomMito = data(1:end,76);
            centerBottomMito = zeros(size(centerBottomMito));
            
            actinwidthmid = data(1:end,85);
            actincentermid = data(1:end,86);
            actinwidthup = data(1:end,91);
            actincenterup = data(1:end,92);
            actinwidthbo = data(1:end,97);
            actincenterbo = data(1:end,98);
            
            if paramPlots
                extraparam = data(1:end,paramCol);
            end
            
            % GET THE MITO WIDTH FITTING DATA AND DECIDE THE WIDTH BASED ON
            % THE THRESHOLD FROM THE FUNCTION INPUTS. OVERWRITE THE MITO
            % WIDTH VECTOR WITH THE NEW WIDTH, AND AS SUCH I CAN KEEP THE
            % REST OF THE CODE AS IT IS. 
            
            % Get the mito mid width fitting data
            widthMidFitData = data(1:end,40:53);
            % Get the mito upper width fitting data
            widthUpperFitData = data(1:end,55:68);
            % Get the mito bottom width fitting data
            widthBottomFitData = data(1:end,70:83);
            
            for i=1:length(areaMito)
                if widthMidFitData(i,2) == 1
                    if widthMidFitData(i,5) > singleTol
                        widthMito(i) = widthMidFitData(i,3);
                        widthMidMito(i) = widthMidFitData(i,3);
                        centerMidMito(i) = widthMidFitData(i,4);
                    else
                        widthMito(i) = 0;
                    end
                elseif widthMidFitData(i,2) == 2
                    if widthMidFitData(i,8) > doubleTol
                        widthMito(i) = widthMidFitData(i,6);
                        widthMidMito(i) = widthMidFitData(i,6);
                        centerMidMito(i) = widthMidFitData(i,7);
                    elseif widthMidFitData(i,5) > doubleTol2
                        widthMito(i) = abs(widthMidFitData(i,9)-widthMidFitData(i,10));
                    else
                        widthMito(i) = 0;
                    end
                else
                    widthMito(i) = 0;
                end
                
                if widthUpperFitData(i,2) == 1
                    if widthUpperFitData(i,5) > singleTol
                        widthUpperMito(i) = widthUpperFitData(i,3);
                        centerUpperMito(i) = widthUpperFitData(i,4);
                    else
                        widthUpperMito(i) = 0;
                    end
                elseif widthUpperFitData(i,2) == 2
                    if widthUpperFitData(i,8) > doubleTol
                        widthUpperMito(i) = widthUpperFitData(i,6);
                        centerUpperMito(i) = widthUpperFitData(i,7);
                    elseif widthUpperFitData(i,5) > doubleTol2
                        widthUpperMito(i) = abs(widthUpperFitData(i,9)-widthUpperFitData(i,10));
                    else
                        widthUpperMito(i) = 0;
                    end
                else
                    widthUpperMito(i) = 0;
                end
                
                if widthBottomFitData(i,2) == 1
                    if widthBottomFitData(i,5) > singleTol
                        widthBottomMito(i) = widthBottomFitData(i,3);
                        centerBottomMito(i) = widthBottomFitData(i,4);
                    else
                        widthBottomMito(i) = 0;
                    end
                elseif widthBottomFitData(i,2) == 2
                    if widthBottomFitData(i,8) > doubleTol
                        widthBottomMito(i) = widthBottomFitData(i,6);
                        centerBottomMito(i) = widthBottomFitData(i,7);
                    elseif widthBottomFitData(i,5) > doubleTol2
                        widthBottomMito(i) = abs(widthBottomFitData(i,9)-widthBottomFitData(i,10));
                    else
                        widthBottomMito(i) = 0;
                    end
                else
                    widthBottomMito(i) = 0;
                end
                
                if widthMito(i) ~= 0 && widthUpperMito(i) ~= 0
                    mitoratioparamum(i) = widthUpperMito(i)/widthMito(i); %Fitted mito upper/mid width ratio
                end
                if widthMito(i) ~= 0 && widthBottomMito(i) ~= 0
                    mitoratioparambm(i) = widthBottomMito(i)/widthMito(i); %Fitted mito bottom/mid width ratio
                end
                if widthUpperMito(i) ~= 0 && widthBottomMito(i) ~= 0
                    mitoratioparamub(i) = min(widthUpperMito(i),widthBottomMito(i))/max(widthUpperMito(i),widthBottomMito(i)); %Fitted mito bottom/upper width ratio
                end
            end
            
            
            for i=1:length(areaMito)
                mitoAreaFiles(i,1,fileNum) = i;
                mitoAreaFiles(i,2,fileNum) = areaMito(i);
                mitoWidthFiles(i,1,fileNum) = i;
                mitoWidthFiles(i,2,fileNum) = widthMito(i);
                mitoCenterFiles(i,1,fileNum) = i;
                mitoCenterFiles(i,2,fileNum) = centerMito(i);
                mitoWidthMidFiles(i,1,fileNum) = i;
                mitoWidthMidFiles(i,2,fileNum) = widthMidMito(i);
                mitoCenterMidFiles(i,1,fileNum) = i;
                mitoCenterMidFiles(i,2,fileNum) = centerMidMito(i);
                mitoWidthUpFiles(i,1,fileNum) = i;
                mitoWidthUpFiles(i,2,fileNum) = widthUpperMito(i);
                mitoCenterUpFiles(i,1,fileNum) = i;
                mitoCenterUpFiles(i,2,fileNum) = centerUpperMito(i);
                mitoWidthBoFiles(i,1,fileNum) = i;
                mitoWidthBoFiles(i,2,fileNum) = widthBottomMito(i);
                mitoCenterBoFiles(i,1,fileNum) = i;
                mitoCenterBoFiles(i,2,fileNum) = centerBottomMito(i);
                % Take the maximum of the ellipsoidal and skeletal length.
                % For small ones, the elliposidal will be longer and is a
                % better match, but for long ones that curve a bit the
                % skeletal fit will hopefully always be better.
                mitoLengthFiles(i,1,fileNum) = i;
                mitoLengthFiles(i,2,fileNum) = max(lengthMito(i),lengthSkelMito(i));
                %mitoLengthSkelFiles(i,1,fileNum) = i;
                %mitoLengthSkelFiles(i,2,fileNum) = lengthSkelMito(i);
%                 mitoNonStripesDistFiles(i,1,fileNum) = i;
%                 mitoNonStripesDistFiles(i,2,fileNum) = nonstripesdistMito(i);
%               % Distance to nearest actin non-stripes binary patch
                mitoNonStripesDistFiles(i,1,fileNum) = i;
                mitoNonStripesDistFiles(i,2,fileNum) = mitopatchdist(i); % Distance to nearest actin patch
                mitoStripesParam(i,1,fileNum) = i;
                mitoStripesParam(i,2,fileNum) = stripesparam(i);
                mitoAISParam(i,1,fileNum) = i;
                mitoAISParam(i,2,fileNum) = aisparam(i);
                mitoSmallParam(i,1,fileNum) = i;
                mitoSmallParam(i,2,fileNum) = smallparam(i);
                mitoBiggerParam(i,1,fileNum) = i;
                mitoBiggerParam(i,2,fileNum) = mitobiggerparam(i);
                mitoRatioParamUMFiles(i,1,fileNum) = i;
                mitoRatioParamUMFiles(i,2,fileNum) = mitoratioparamum(i);
                mitoRatioParamBMFiles(i,1,fileNum) = i;
                mitoRatioParamBMFiles(i,2,fileNum) = mitoratioparambm(i);
                mitoRatioParamUBFiles(i,1,fileNum) = i;
                mitoRatioParamUBFiles(i,2,fileNum) = mitoratioparamub(i);
                filamentWidthFiles(i,1,fileNum) = i;
                filamentWidthFiles(i,2,fileNum) = actinwidth(i);
                filamentCenterFiles(i,1,fileNum) = i;
                filamentCenterFiles(i,2,fileNum) = actincenter(i);
                mitoSmallAreaParamFiles(i,1,fileNum) = i;
                mitoSmallAreaParamFiles(i,2,fileNum) = smallareaparam(i);
                filamentDoublePeakParamFiles(i,1,fileNum) = i;
                filamentDoublePeakParamFiles(i,2,fileNum) = actindoublepeakparam(i);
                filamentSinglePeakParamFiles(i,1,fileNum) = i;
                filamentSinglePeakParamFiles(i,2,fileNum) = actinsinglepeakparam(i);
                mitoActinDistFiles(i,1,fileNum) = i;
                mitoActinDistFiles(i,2,fileNum) = mitoactindist(i);
                mitoActinDistNegParamFiles(i,1,fileNum) = i;
                mitoActinDistNegParamFiles(i,2,fileNum) = mitoactindistnegparam(i);
                mitoActinDistD2Files(i,1,fileNum) = i;
                mitoActinDistD2Files(i,2,fileNum) = mitoactindistD2(i);
                mitoActinDistD2CloseParamFiles(i,1,fileNum) = i;
                mitoActinDistD2CloseParamFiles(i,2,fileNum) = mitoactindistD2closeparam(i);
                mitoPatchOverlapFiles(i,1,fileNum) = i;
                mitoPatchOverlapFiles(i,2,fileNum) = mitopatchoverlap(i);
                filamentWidthMidFiles(i,1,fileNum) = i;
                filamentWidthMidFiles(i,2,fileNum) = actinwidthmid(i);
                filamentCenterMidFiles(i,1,fileNum) = i;
                filamentCenterMidFiles(i,2,fileNum) = actincentermid(i);
                filamentWidthUpFiles(i,1,fileNum) = i;
                filamentWidthUpFiles(i,2,fileNum) = actinwidthup(i);
                filamentCenterUpFiles(i,1,fileNum) = i;
                filamentCenterUpFiles(i,2,fileNum) = actincenterup(i);
                filamentWidthBoFiles(i,1,fileNum) = i;
                filamentWidthBoFiles(i,2,fileNum) = actinwidthbo(i);
                filamentCenterBoFiles(i,1,fileNum) = i;
                filamentCenterBoFiles(i,2,fileNum) = actincenterbo(i);
                if paramPlots 
                    extraParamFiles(i,1,fileNum) = i;
                    extraParamFiles(i,2,fileNum) = extraparam(i);
                end
            end
        catch err
            disp(strcat(num2str(fileNum),': No image with this number or a file reading error.'))
        end 
    end

    mitoWidthS = [];
    mitoLengthS = [];
    mitoAreaS = [];
    mitoCenterS = [];
    mitoNonStripesDistS = [];
    filamentWidthS = [];
    filamentCenterS = [];
    mitoDisplacementEdgeS = [];
    mitoWidthRatioUMS = [];
    mitoWidthRatioBMS = [];
    mitoWidthRatioUBS = [];
    mitoActinDistS = [];
    mitoActinDistD2S = [];
    mitoPatchOverlapS = [];

    mitoWidthNS = [];
    mitoLengthNS = [];
    mitoAreaNS = [];
    mitoCenterNS = [];
    mitoNonStripesDistNS = [];
    filamentWidthNS = [];
    filamentCenterNS = [];
    mitoDisplacementEdgeNS = [];
    mitoWidthRatioUMNS = [];
    mitoWidthRatioBMNS = [];
    mitoWidthRatioUBNS = [];
    mitoActinDistNS = [];
    mitoActinDistD2NS = [];
    mitoPatchOverlapNS = [];

    mitoWidthP = [];
    mitoLengthP = [];
    mitoAreaP = [];
    mitoCenterP = [];
    mitoNonStripesDistP = [];
    filamentWidthP = [];
    filamentCenterP = [];
    mitoDisplacementEdgeP = [];
    mitoWidthRatioUMP = [];
    mitoWidthRatioBMP = [];
    mitoWidthRatioUBP = [];
    mitoActinDistP = [];
    mitoActinDistD2P = [];
    mitoPatchOverlapP = [];

    mitoWidthNP = [];
    mitoLengthNP = [];
    mitoAreaNP = [];
    mitoCenterNP = [];
    mitoNonStripesDistNP = [];
    filamentWidthNP = [];
    filamentCenterNP = [];
    mitoDisplacementEdgeNP = [];
    mitoWidthRatioUMNP = [];
    mitoWidthRatioBMNP = [];
    mitoWidthRatioUBNP = [];
    mitoActinDistNP = [];
    mitoActinDistD2NP = [];
    mitoPatchOverlapNP = [];

    mitoThreePeaksNonStripesDist = [];
    noOfThreePeaksInNonStripes = 0;
    noOfThreePeaksInStripes = 0;

    %%% DISPLAY THE MITOCHONDRIA NUMBER OF ALL THE MITOCHONDRIA WITH A WIDTH
    %%% RATIO SMALLER THAN THE SPECIFIC NUMBER. 

    disp(['All mitochondria with at least one width ratio below ' num2str(mitoWidthRatioThresholdLow)])
    for fileNum=fileNumbers
        for i=1:mitosPerFile
            if mitoRatioParamUMFiles(i,2,fileNum) ~= 0 && mitoRatioParamUMFiles(i,2,fileNum) < mitoWidthRatioThresholdLow
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i), ', U/M'])
            end
            if mitoRatioParamBMFiles(i,2,fileNum) ~= 0 && mitoRatioParamBMFiles(i,2,fileNum) < mitoWidthRatioThresholdLow
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i), ', B/M'])
            end
            if mitoRatioParamUBFiles(i,2,fileNum) ~= 0 && mitoRatioParamUBFiles(i,2,fileNum) < mitoWidthRatioThresholdLow
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i), ', U/B'])
            end
        end
    end

    disp(['All mitochondria with at least one width ratio above ' num2str(mitoWidthRatioThresholdHigh)])
    for fileNum=fileNumbers
        for i=1:mitosPerFile
            if mitoRatioParamUMFiles(i,2,fileNum) ~= 0 && mitoRatioParamUMFiles(i,2,fileNum) > mitoWidthRatioThresholdHigh
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i), ', U/M'])
            end
            if mitoRatioParamBMFiles(i,2,fileNum) ~= 0 && mitoRatioParamBMFiles(i,2,fileNum) > mitoWidthRatioThresholdHigh
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i), ', B/M'])
            end
        end
    end

    if paramPlots
        disp(['All mitochondria interesting for the given extra parameter: ' paramName])
    end

    for fileNum=fileNumbers
        noProfiles = size(mitoThreePeaksNonStripesDistFiles);
        for i=1:noProfiles(1)
            if mitoThreePeaksWidthOneFiles(i,2,fileNum) ~= 0
                mitoThreePeaksNonStripesDist = [mitoThreePeaksNonStripesDist; mitoThreePeaksNonStripesDistFiles(i,2,fileNum)];
                if mitoThreePeaksNonStripesDistFiles(i,2,fileNum) == 0
                    noOfThreePeaksInNonStripes = noOfThreePeaksInNonStripes+1;
                else
                    noOfThreePeaksInStripes = noOfThreePeaksInStripes+1;
                end
            end
        end

        % OLD MITO DISPLACEMENT CALCULATION, ONE-FORMULA. WILL CACLULATE
        % ERRONEOUSLY FOR A LOT OF MITOCHONDRIA. ACTINDIST IS THE NEW
        % FORMULA. 
        for i=1:mitosPerFile
            mitoDisplacementEdge = NaN(1,3);
            if filamentWidthMidFiles(i,2,fileNum) ~= 0 && mitoWidthMidFiles(i,2,fileNum) ~= 0
                mitoDisplacementEdge(1) = filamentWidthMidFiles(i,2,fileNum)/2-abs((filamentCenterMidFiles(i,2,fileNum)-actinLineProfLength/2)-(mitoCenterMidFiles(i,2,fileNum)-mitoLineProfLength/2))+mitoWidthMidFiles(i,2,fileNum)/2;
            end
            if filamentWidthUpFiles(i,2,fileNum) ~= 0 && mitoWidthUpFiles(i,2,fileNum) ~= 0
                mitoDisplacementEdge(2) = filamentWidthUpFiles(i,2,fileNum)/2-abs((filamentCenterUpFiles(i,2,fileNum)-actinLineProfLength/2)-(mitoCenterUpFiles(i,2,fileNum)-mitoLineProfLength/2))+mitoWidthUpFiles(i,2,fileNum)/2;
            end
            if filamentWidthBoFiles(i,2,fileNum) ~= 0 && mitoWidthBoFiles(i,2,fileNum) ~= 0
                mitoDisplacementEdge(3) = filamentWidthBoFiles(i,2,fileNum)/2-abs((filamentCenterBoFiles(i,2,fileNum)-actinLineProfLength/2)-(mitoCenterBoFiles(i,2,fileNum)-mitoLineProfLength/2))+mitoWidthBoFiles(i,2,fileNum)/2;
            end
            
            if mitoAreaFiles(i,2,fileNum) ~= 0
                if mitoStripesParam(i,2,fileNum) == 1
                    mitoPatchOverlapS = [mitoPatchOverlapS; mitoPatchOverlapFiles(i,2,fileNum)];
%                     mitoAreaS = [mitoAreaS; mitoAreaFiles(i,2,fileNum)];
%                     mitoLengthS = [mitoLengthS; mitoLengthFiles(i,2,fileNum)];
%                     mitoCenterS = [mitoCenterS; mitoCenterFiles(i,2,fileNum)];
                    if mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoWidthS = [mitoWidthS; mitoWidthFiles(i,2,fileNum)];
                        mitoAreaS = [mitoAreaS; mitoAreaFiles(i,2,fileNum)];
                        mitoLengthS = [mitoLengthS; mitoLengthFiles(i,2,fileNum)];
                        mitoCenterS = [mitoCenterS; mitoCenterFiles(i,2,fileNum)];
                        mitoNonStripesDistS = [mitoNonStripesDistS; mitoNonStripesDistFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0
                        filamentWidthS = [filamentWidthS; filamentWidthFiles(i,2,fileNum)];
                        filamentCenterS = [filamentCenterS; filamentCenterFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamUMFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioUMS = [mitoWidthRatioUMS; mitoRatioParamUMFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamBMFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioBMS = [mitoWidthRatioBMS; mitoRatioParamBMFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamUBFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioUBS = [mitoWidthRatioUBS; mitoRatioParamUBFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0 && mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoActinDistS = [mitoActinDistS; mitoActinDistFiles(i,2,fileNum)];
                        mitoActinDistD2S = [mitoActinDistD2S; mitoActinDistD2Files(i,2,fileNum)];
                        mitoDisplacementEdgeS = [mitoDisplacementEdgeS; min(mitoDisplacementEdge)];
                    end
                elseif mitoStripesParam(i,2,fileNum) == 0
                    mitoPatchOverlapNS = [mitoPatchOverlapNS; mitoPatchOverlapFiles(i,2,fileNum)];
%                     mitoAreaNS = [mitoAreaNS; mitoAreaFiles(i,2,fileNum)];
%                     mitoLengthNS = [mitoLengthNS; mitoLengthFiles(i,2,fileNum)];
%                     mitoCenterNS = [mitoCenterNS; mitoCenterFiles(i,2,fileNum)];
                    if mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoWidthNS = [mitoWidthNS; mitoWidthFiles(i,2,fileNum)];
                        mitoAreaNS = [mitoAreaNS; mitoAreaFiles(i,2,fileNum)];
                        mitoLengthNS = [mitoLengthNS; mitoLengthFiles(i,2,fileNum)];
                        mitoCenterNS = [mitoCenterNS; mitoCenterFiles(i,2,fileNum)];
                        mitoNonStripesDistNS = [mitoNonStripesDistNS; mitoNonStripesDistFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0
                        filamentWidthNS = [filamentWidthNS; filamentWidthFiles(i,2,fileNum)];
                        filamentCenterNS = [filamentCenterNS; filamentCenterFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamUMFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioUMNS = [mitoWidthRatioUMNS; mitoRatioParamUMFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamBMFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioBMNS = [mitoWidthRatioBMNS; mitoRatioParamBMFiles(i,2,fileNum)];
                    end
                    if mitoRatioParamUBFiles(i,2,fileNum) ~= 0
                        mitoWidthRatioUBNS = [mitoWidthRatioUBNS; mitoRatioParamUBFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0 && mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoActinDistNS = [mitoActinDistNS; mitoActinDistFiles(i,2,fileNum)];
                        mitoActinDistD2NS = [mitoActinDistD2NS; mitoActinDistD2Files(i,2,fileNum)];
                        mitoDisplacementEdgeNS = [mitoDisplacementEdgeNS; min(mitoDisplacementEdge)];
                    end
                end
                if extraParamFiles(i,2,fileNum) == 1
                    %disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i)])
                    mitoPatchOverlapP = [mitoPatchOverlapP; mitoPatchOverlapFiles(i,2,fileNum)];
%                     mitoAreaP = [mitoAreaP; mitoAreaFiles(i,2,fileNum)];
%                     mitoLengthP = [mitoLengthP; mitoLengthFiles(i,2,fileNum)];
%                     mitoCenterP = [mitoCenterP; mitoCenterFiles(i,2,fileNum)];
                    if mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoWidthP = [mitoWidthP; mitoWidthFiles(i,2,fileNum)];
                        mitoAreaP = [mitoAreaP; mitoAreaFiles(i,2,fileNum)];
                        mitoLengthP = [mitoLengthP; mitoLengthFiles(i,2,fileNum)];
                        mitoCenterP = [mitoCenterP; mitoCenterFiles(i,2,fileNum)];
                        mitoNonStripesDistP = [mitoNonStripesDistP; mitoNonStripesDistFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0
                        filamentWidthP = [filamentWidthP; filamentWidthFiles(i,2,fileNum)];
                        filamentCenterP = [filamentCenterP; filamentCenterFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0 && mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoActinDistP = [mitoActinDistP; mitoActinDistFiles(i,2,fileNum)];
                        mitoActinDistD2P = [mitoActinDistD2P; mitoActinDistD2Files(i,2,fileNum)];
                        mitoDisplacementEdgeP = [mitoDisplacementEdgeP; min(mitoDisplacementEdge)];
                    end
                elseif extraParamFiles(i,2,fileNum) == 0
                    mitoPatchOverlapNP = [mitoPatchOverlapNP; mitoPatchOverlapFiles(i,2,fileNum)];
%                     mitoAreaNP = [mitoAreaNP; mitoAreaFiles(i,2,fileNum)];
%                     mitoLengthNP = [mitoLengthNP; mitoLengthFiles(i,2,fileNum)];
%                     mitoCenterNP = [mitoCenterNP; mitoCenterFiles(i,2,fileNum)];
                    if mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoWidthNP = [mitoWidthNP; mitoWidthFiles(i,2,fileNum)];
                        mitoAreaNP = [mitoAreaNP; mitoAreaFiles(i,2,fileNum)];
                        mitoLengthNP = [mitoLengthNP; mitoLengthFiles(i,2,fileNum)];
                        mitoCenterNP = [mitoCenterNP; mitoCenterFiles(i,2,fileNum)];
                        mitoNonStripesDistNP = [mitoNonStripesDistNP; mitoNonStripesDistFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0
                        filamentWidthNP = [filamentWidthNP; filamentWidthFiles(i,2,fileNum)];
                        filamentCenterNP = [filamentCenterNP; filamentCenterFiles(i,2,fileNum)];
                    end
                    if filamentWidthFiles(i,2,fileNum) ~= 0 && mitoWidthFiles(i,2,fileNum) ~= 0
                        mitoActinDistNP = [mitoActinDistNP; mitoActinDistFiles(i,2,fileNum)];
                        mitoActinDistD2NP = [mitoActinDistD2NP; mitoActinDistD2Files(i,2,fileNum)];
                        mitoDisplacementEdgeNP = [mitoDisplacementEdgeNP; min(mitoDisplacementEdge)];
                    end
                end
            end
        end 
    end

    [nonzeroS,~,filamentWidthSNZ]=find(filamentWidthS);
    [nonzeroNS,~,filamentWidthNSNZ]=find(filamentWidthNS);
    [nonzeroP,~,filamentWidthPNZ]=find(filamentWidthP);
    [nonzeroNP,~,filamentWidthNPNZ]=find(filamentWidthNP);
    
%     assignin('base','mitowidthstedGT14',[mitoWidthS;mitoWidthNS]);
%     assignin('base','mitolengthstedGT14',[mitoLengthS;mitoLengthNS]);
%     assignin('base','mitowidthstedS',mitoWidthS);
%     assignin('base','mitowidthstedNS',mitoWidthNS);
%     assignin('base','mitolengthstedS',mitoLengthS);
%     assignin('base','mitolengthstedNS',mitoLengthNS);
%     assignin('base','mitowidthsted',[mitoWidthS;mitoWidthNS]);
%     assignin('base','mitolengthsted',[mitoLengthS;mitoLengthNS]);
%     assignin('base','mitoareasted',[mitoAreaS;mitoAreaNS]);
%     assignin('base','mitowidthratio2',[mitoWidthRatioUMS;mitoWidthRatioUMNS;mitoWidthRatioBMS;mitoWidthRatioBMNS]);
%     assignin('base','mitopatchoverlap',[mitoPatchOverlapS;mitoPatchOverlapNS]);
%     assignin('base','mitopatchdiststedS',mitoNonStripesDistS);
%     assignin('base','mitopatchdiststedNS',mitoNonStripesDistNS);
%     assignin('base','actinwidthstedS2',filamentWidthS);
%     assignin('base','actinwidthstedNS2',filamentWidthNS);
%      assignin('base','actinmitodiststedS',mitoActinDistS);
%      assignin('base','actinmitodiststedNS',mitoActinDistNS);
%      assignin('base','mitowidthBigFOVDIV17Ax2',mitoWidthP);
%      assignin('base','mitowidthBigFOVDIV17De2',mitoWidthNP);
%      assignin('base','mitolengthOnlyFitBigFOVDIV17Ax2',mitoLengthP);
%      assignin('base','mitolengthOnlyFitBigFOVDIV17De2',mitoLengthNP);
%      assignin('base','mitoareaBigFOVDIV17Ax2',mitoAreaP);
%      assignin('base','mitoareaBigFOVDIV17De2',mitoAreaNP);
%      assignin('base','actinwidthBigFOVDIV17Ax2',filamentWidthP);
%      assignin('base','actinwidthBigFOVDIV17De2',filamentWidthNP);
%      assignin('base','actinwidthBigFOVDIV3Ax',filamentWidthP);
%      assignin('base','actinwidthBigFOVDIV3De',filamentWidthNP);
%      assignin('base','mitowidthBigFOVDIV3Ax',mitoWidthP);
%      assignin('base','mitowidthBigFOVDIV3De',mitoWidthNP);
%      assignin('base','mitolengthOnlyFitBigFOVDIV3Ax',mitoLengthP);
%      assignin('base','mitolengthOnlyFitBigFOVDIV3De',mitoLengthNP);
%      assignin('base','mitoareaBigFOVDIV3Ax',mitoAreaP);
%      assignin('base','mitoareaBigFOVDIV3De',mitoAreaNP);
%      assignin('base','actinwidthAISSomaAx',filamentWidthP);
%      assignin('base','actinwidthAISSomaDe',filamentWidthNP);
%      assignin('base','mitowidthAISSomaAx',mitoWidthP);
%      assignin('base','mitowidthAISSomaDe',mitoWidthNP);
%      assignin('base','mitolengthOnlyFitAISSomaAx',mitoLengthP);
%      assignin('base','mitolengthOnlyFitAISSomaDe',mitoLengthNP);
%      assignin('base','mitoareaAISSomaAx',mitoAreaP);
%      assignin('base','mitoareaAISSomaDe',mitoAreaNP);
%      assignin('base','actinmitodistD2P',mitoActinDistD2P);
%      assignin('base','actinmitodistD2NP',mitoActinDistD2NP);
%     assignin('base','mitowidthstedfixedS',mitoWidthS);
%     assignin('base','mitowidthstedfixedNS',mitoWidthNS);
%     assignin('base','mitolengthstedfixedS',mitoLengthS);
%     assignin('base','mitolengthstedfixedNS',mitoLengthNS);
%     assignin('base','mitoareastedfixedS',mitoAreaS);
%     assignin('base','mitoareastedfixedNS',mitoAreaNS);
%     assignin('base','mitowidthratiostedUBS',mitoWidthRatioUBS);
%     assignin('base','mitowidthratiostedUBNS',mitoWidthRatioUBNS);
%     assignin('base','mitowidthratiostedUMS',mitoWidthRatioUMS);
%     assignin('base','mitowidthratiostedUMNS',mitoWidthRatioUMNS);
%     assignin('base','mitowidthratiostedBMS',mitoWidthRatioBMS);
%     assignin('base','mitowidthratiostedBMNS',mitoWidthRatioBMNS);
%     assignin('base','mitowidthratiostedall',[mitoWidthRatioUMS;mitoWidthRatioUMNS;mitoWidthRatioBMS;mitoWidthRatioBMNS]);
%     assignin('base','mitoareastedS',mitoAreaS);
%     assignin('base','mitoareastedNS',mitoAreaNS);

%%% AGES DIVGT9 & DIVLT9 & DIVLT4
%     assignin('base','actinwidthGT9',[filamentWidthS;filamentWidthNS]);
%     assignin('base','mitowidthGT9',[mitoWidthS;mitoWidthNS]);
%     assignin('base','mitolengthGT9',[mitoLengthS;mitoLengthNS]);
%     assignin('base','mitoareaGT9',[mitoAreaS;mitoAreaNS]);
%     assignin('base','mitowidthratioGT9',[mitoWidthRatioUMS;mitoWidthRatioUMNS;mitoWidthRatioBMS;mitoWidthRatioBMNS]);
%     assignin('base','actinmitodistD2GT9',[mitoActinDistD2S;mitoActinDistD2NS]);
%     assignin('base','mitopatchoverlapGT9',[mitoPatchOverlapS;mitoPatchOverlapNS]);
%     assignin('base','mitopatchdiststedGT9',[mitoNonStripesDistS;mitoNonStripesDistNS]);
%     
%     assignin('base','actinwidthGT9S',filamentWidthS);
%     assignin('base','mitowidthGT9S',mitoWidthS);
%     assignin('base','mitolengthGT9S',mitoLengthS);
%     assignin('base','mitoareaGT9S',mitoAreaS);
%     assignin('base','mitowidthratioGT9S',[mitoWidthRatioUMS;mitoWidthRatioBMS]);
%     assignin('base','actinmitodistD2GT9S',mitoActinDistD2S);
%     assignin('base','mitopatchoverlapGT9S',mitoPatchOverlapS);
%     assignin('base','mitopatchdiststedGT9S',mitoNonStripesDistS);
%     
%     assignin('base','actinwidthGT9NS',filamentWidthNS);
%     assignin('base','mitowidthGT9NS',mitoWidthNS);
%     assignin('base','mitolengthGT9NS',mitoLengthNS);
%     assignin('base','mitoareaGT9NS',mitoAreaNS);
%     assignin('base','mitowidthratioGT9NS',[mitoWidthRatioUMNS;mitoWidthRatioBMNS]);
%     assignin('base','actinmitodistD2GT9NS',mitoActinDistD2NS);
%     assignin('base','mitopatchoverlapGT9NS',mitoPatchOverlapNS);
%     assignin('base','mitopatchdiststedGT9NS',mitoNonStripesDistNS);

%     assignin('base','actinwidthLT4',[filamentWidthS;filamentWidthNS]);
%     assignin('base','mitowidthLT4',[mitoWidthS;mitoWidthNS]);
%     assignin('base','mitolengthLT4',[mitoLengthS;mitoLengthNS]);
%     assignin('base','mitoareaLT4',[mitoAreaS;mitoAreaNS]);
%     assignin('base','mitowidthratioLT4',[mitoWidthRatioUMS;mitoWidthRatioUMNS;mitoWidthRatioBMS;mitoWidthRatioBMNS]);
%     assignin('base','actinmitodistD2LT4',[mitoActinDistD2S;mitoActinDistD2NS]);
%     assignin('base','mitopatchoverlapLT4',[mitoPatchOverlapS;mitoPatchOverlapNS]);
%     assignin('base','mitopatchdiststedLT4',[mitoNonStripesDistS;mitoNonStripesDistNS]);
%     
%     assignin('base','actinwidthLT4S',filamentWidthS);
%     assignin('base','mitowidthLT4S',mitoWidthS);
%     assignin('base','mitolengthLT4S',mitoLengthS);
%     assignin('base','mitoareaLT4S',mitoAreaS);
%     assignin('base','mitowidthratioLT4S',[mitoWidthRatioUMS;mitoWidthRatioBMS]);
%     assignin('base','actinmitodistD2LT4S',mitoActinDistD2S);
%     assignin('base','mitopatchoverlapLT4S',mitoPatchOverlapS);
%     assignin('base','mitopatchdiststedLT4S',mitoNonStripesDistS);
%     
%     assignin('base','actinwidthLT4NS',filamentWidthNS);
%     assignin('base','mitowidthLT4NS',mitoWidthNS);
%     assignin('base','mitolengthLT4NS',mitoLengthNS);
%     assignin('base','mitoareaLT4NS',mitoAreaNS);
%     assignin('base','mitowidthratioLT4NS',[mitoWidthRatioUMNS;mitoWidthRatioBMNS]);
%     assignin('base','actinmitodistD2LT4NS',mitoActinDistD2NS);
%     assignin('base','mitopatchoverlapLT4NS',mitoPatchOverlapNS);
%     assignin('base','mitopatchdiststedLT4NS',mitoNonStripesDistNS);

    
    assignin('base','actinwidthGT4',[filamentWidthS;filamentWidthNS]);
    assignin('base','mitowidthGT4',[mitoWidthS;mitoWidthNS]);
    widsave = [mitoWidthS;mitoWidthNS];
    assignin('base','mitolengthGT4',[mitoLengthS;mitoLengthNS]);
    lensave = [mitoLengthS;mitoLengthNS];
    assignin('base','mitoareaGT4',[mitoAreaS;mitoAreaNS]);
    aresave = [mitoAreaS;mitoAreaNS];
    assignin('base','mitowidthratioGT4',[mitoWidthRatioUMS;mitoWidthRatioUMNS;mitoWidthRatioBMS;mitoWidthRatioBMNS]);
    assignin('base','actinmitodistD2GT4',[mitoActinDistD2S;mitoActinDistD2NS]);
    assignin('base','mitopatchoverlapGT4',[mitoPatchOverlapS;mitoPatchOverlapNS]);
    assignin('base','mitopatchdiststedGT4',[mitoNonStripesDistS;mitoNonStripesDistNS]);
    
    assignin('base','actinwidthGT4S',filamentWidthS);
    assignin('base','mitowidthGT4S',mitoWidthS);
    assignin('base','mitolengthGT4S',mitoLengthS);
    assignin('base','mitoareaGT4S',mitoAreaS);
    assignin('base','mitowidthratioGT4S',[mitoWidthRatioUMS;mitoWidthRatioBMS]);
    assignin('base','actinmitodistD2GT4S',mitoActinDistD2S);
    assignin('base','mitopatchoverlapGT4S',mitoPatchOverlapS);
    assignin('base','mitopatchdiststedGT4S',mitoNonStripesDistS);
    
    assignin('base','actinwidthGT4NS',filamentWidthNS);
    assignin('base','mitowidthGT4NS',mitoWidthNS);
    assignin('base','mitolengthGT4NS',mitoLengthNS);
    assignin('base','mitoareaGT4NS',mitoAreaNS);
    assignin('base','mitowidthratioGT4NS',[mitoWidthRatioUMNS;mitoWidthRatioBMNS]);
    assignin('base','actinmitodistD2GT4NS',mitoActinDistD2NS);
    assignin('base','mitopatchoverlapGT4NS',mitoPatchOverlapNS);
    assignin('base','mitopatchdiststedGT4NS',mitoNonStripesDistNS);
    
    patchdistPsave = mitoNonStripesDistP;
    patchdistNPsave = mitoNonStripesDistNP;
    patchoverlapPsave = mitoPatchOverlapP;
    patchoverlapNPsave = mitoPatchOverlapNP;
    
%     assignin('base','mitowidthstedfixed',mitoWidthS);
%     assignin('base','mitolengthstedfixed',mitoLengthS);
%     assignin('base','mitoareastedfixed',mitoAreaS);
%     assignin('base','mitowidthratiostedfixed',[mitoWidthRatioUMS;mitoWidthRatioBMS]);
%     assignin('base','mitopatchoverlapstedfixed',mitoPatchOverlapS);
%     assignin('base','mitopatchdiststedstedfixed',mitoNonStripesDistS);
%     assignin('base','actinwidthstedstedfixed',filamentWidthS);
%     assignin('base','actinmitodiststedfixed',mitoActinDistD2S);
    save('tempVar.mat','widsave','lensave','aresave','patchdistPsave','patchdistNPsave','patchoverlapPsave','patchoverlapNPsave');
end