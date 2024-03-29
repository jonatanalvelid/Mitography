%%% 
%
% Get all the results from the mitography, into variables for morphology
% and boolean variables.
% Version 201203 - OXPHOS
%
%%%
clear

% Add function folder to filepath, so that those functions can be read.
functionFolder = fileparts(which('findFunctionFolders.m'));
addpath(genpath(functionFolder));

masterFolderPath = strcat(uigetdir('E:\PhD\data_analysis\antimycin'),'\');
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

areathresh_wEll = 0.1;
mitosPerFile = 1000;

lastFileNumber = 50;

fileNumbersCells4 = [1 3;4 5;6 7;8 10;11 13;14 16;17 19;20 22;23 25;26 27;28 30;31 33;34 36;37 39;40 41];
fileNumbersCells5 = [1 2;3 6;7 9;10 13;14 16;17 19;20 22;23 25;26 27;28 30;31 33;34 35;36 37;38 39;41 42;43 44;45 47;48 49];

filenameAnalysis = '_MitoAnalysisFull.txt';

%%% TAKE CARE OF ALL MITOCHONDRIA DATA, DIVIDE FOR AIS AND STRIPES, AND
%%% THREE PEAKS

expnumber = 4;
ctaa = 'ct';
fileNumbersCells = fileNumbersCells4;
nocells = size(fileNumbersCells);
nocells = nocells(1);
dataset = split(masterFolderPath,'\');
dataset = char(dataset(5));

figure()
subplot(ceil(nocells/5),5,1)

for cellno = 1:nocells
    fileNumbers = fileNumbersCells(cellno,1):fileNumbersCells(cellno,2);
    
    mitoWidthFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoWidthEllFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoLengthFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitoAreaFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitosoxParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    mitosoxvalFiles = zeros(mitosPerFile,2,lastFileNumber);
    ompvalFiles = zeros(mitosPerFile,2,lastFileNumber);
    somaParamFiles = zeros(mitosPerFile,2,lastFileNumber);
    borderParamFiles = zeros(mitosPerFile,2,lastFileNumber);

    for fileNum = fileNumbers
        if fileNum < 10
            filename = strcat('Image_00',int2str(fileNum),filenameAnalysis);
        else
            filename = strcat('Image_0',int2str(fileNum),filenameAnalysis);
        end
        filepath = strcat(masterFolderPath,filename);

        try
            data = dlmread(filepath);
            %disp(filename);
            areaMito = data(1:end,4);
            lengthMito = data(1:end,5); %Ellipsoidal fit mitochondria length (major axis)
            widthEllMito = data(1:end,6); %Ellipsoidal fit mitochondria width (minor axis)
            lengthSkelMito = data(1:end,7); %Skeleton mitochondria length (skeleton part closest to the mitochondria centroid)
            mitosoxval = data(1:end,111);
            mitosoxparam = data(1:end,112);
            ompval = data(1:end,115);

            somaparam = data(1:end,109);
            borderparam = data(1:end,110);

            for i=1:length(areaMito)
                mitoAreaFiles(i,1,fileNum) = i;
                mitoAreaFiles(i,2,fileNum) = areaMito(i);
                mitoWidthEllFiles(i,1,fileNum) = i;
                mitoWidthEllFiles(i,2,fileNum) = widthEllMito(i);
                % Take the maximum of the ellipsoidal and skeletal length.
                % For small ones, the elliposidal will be longer and is a
                % better match, but for long ones that curve a bit or branches
                % the skeletal fit will hopefully always be better. However,
                % sometimes for small ones, when no skeleton exists, it will
                % pick the close-by skeleton, which might be a big one, or
                % entirely the wrong length. Thus always pick the ellipsoidal
                % length for these. 
                mitoLengthFiles(i,1,fileNum) = i;
                if areaMito(i) < 0.2
                    mitoLengthFiles(i,2,fileNum) = lengthMito(i);
                else
                    mitoLengthFiles(i,2,fileNum) = max(lengthMito(i),lengthSkelMito(i));
                end

                mitosoxParamFiles(i,1,fileNum) = i;
                mitosoxParamFiles(i,2,fileNum) = mitosoxparam(i);
                mitosoxvalFiles(i,1,fileNum) = i;
                mitosoxvalFiles(i,2,fileNum) = mitosoxval(i);
                ompvalFiles(i,1,fileNum) = i;
                ompvalFiles(i,2,fileNum) = ompval(i);

                somaParamFiles(i,1,fileNum) = i;
                somaParamFiles(i,2,fileNum) = somaparam(i);
                borderParamFiles(i,1,fileNum) = i;
                borderParamFiles(i,2,fileNum) = borderparam(i);

            end
        catch err
            disp(strcat(num2str(fileNum),': No image with this number or a file reading error.'))
        end  
    end

    mitoWidth = [];
    mitoLength = [];
    mitoAR = [];
    mitoArea = [];
    mitoOXPHOS = [];
    mitoOMP = [];
    mitoOXPHOSparam = [];
    mitodoublepeakparam = [];

    mitoWidthP = [];
    mitoLengthP = [];
    mitoARP = [];
    mitoAreaP = [];
    % mitoTMREP = [];
    % mitoOMPP = [];
    % mitoTMREparamP = [];
    mitodoublepeakparamP = [];

    mitoWidthNP = [];
    mitoLengthNP = [];
    mitoARNP = [];
    mitoAreaNP = [];
    % mitoTMRENP = [];
    % mitoOMPNP = [];
    % mitoTMREparamNP = [];
    mitodoublepeakparamNP = [];


    for fileNum=fileNumbers
        for i=1:mitosPerFile
            allcheck = somaParamFiles(i,2,fileNum) | borderParamFiles(i,2,fileNum);
            if mitoWidthEllFiles(i,2,fileNum) ~= 0 && ~allcheck
                % Calculate AR as w_ell/l_ell if the area is small enough
                % (A<0.2 µm^2), while instead use w_fit/l_ell if the
                % mitochondria is bigger. The fitted width will always be the
                % more accurate width, but since we don't have a completely 
                % accurate (fitted) length for the small mitos, the AR will be 
                % more accurate by using the ellipsoidal width and length for 
                % the smaller mitos. 
                ARtemp = mitoWidthEllFiles(i,2,fileNum)/mitoLengthFiles(i,2,fileNum);
                if ARtemp > 1
                    ARtemp = 1/ARtemp;
                end
                mitoWidth = [mitoWidth; mitoWidthEllFiles(i,2,fileNum)];
                mitoArea = [mitoArea; mitoAreaFiles(i,2,fileNum)];
                mitoLength = [mitoLength; mitoLengthFiles(i,2,fileNum)];
                mitoAR = [mitoAR; ARtemp];
                mitoOXPHOS = [mitoOXPHOS; mitosoxvalFiles(i,2,fileNum)];
                mitoOMP = [mitoOMP; ompvalFiles(i,2,fileNum)];
                mitoOXPHOSparam = [mitoOXPHOSparam; mitosoxParamFiles(i,2,fileNum)];
            end
        end 
    end
    
    %disp(size(mitoArea))
    %disp(mitoArea)

    %{
    %%% DISPLAY THE MITOCHONDRIA NUMBER OF SPECIFIC MITOCHONDRIA WITH A CERTAIN
    %%% PARAMETER ABOVE/BELOW A CERTAIN THRESHOLD
    mitoAreaThresh = 0.086;
    mitoARThresh = 0.5;

    disp(['All stick-like small mitos with TMRE+: '])
    for fileNum=fileNumbers
        for i=1:mitosPerFile
            ARtemp = mitoWidthFiles(i,2,fileNum)/mitoLengthFiles(i,2,fileNum);
            if ARtemp > 1
                ARtemp = 1/ARtemp;
            end
            if tmreParamFiles(i,2,fileNum) == 0 && mitoAreaFiles(i,2,fileNum) < mitoAreaThresh && ARtemp < mitoARThresh
                disp(['Image: ' int2str(fileNum) ', mito: ' int2str(i)])
            end
        end
    end
    %}

    subplot(ceil(nocells/5),5,cellno)
    hold on
    %histogram(mitoArea(mitoArea<0.086),'BinWidth',0.01)
    histogram(mitoArea(mitoArea<0.086 & mitoOXPHOSparam==0),'FaceColor','b','BinWidth',0.01)
    histogram(mitoArea(mitoArea<0.086 & mitoOXPHOSparam==1),'FaceColor','r','BinWidth',0.01)
    title(sprintf('%s - cellno: %i',dataset,cellno))
    legend({'OXPHOS-','OXPHOS+'})
    xlim([0 0.1])
    
    mitoinfo = struct();
    if ~isempty(mitoArea)
        mitoinfo.exp = expnumber;
        mitoinfo.ctaa = ctaa;
        mitoinfo.cell = cellno;
        mitoinfo.area = mitoArea;
        mitoinfo.ar = mitoAR;
        mitoinfo.length = mitoLength;
        mitoinfo.width = mitoWidth;
        mitoinfo.omp = mitoOMP;
        mitoinfo.oxphos = mitoOXPHOS;
        mitoinfo.oxphosparam = mitoOXPHOSparam;
        savename = sprintf('oxphos_e%i_%s_c%i',expnumber,ctaa,cellno);
        save(savename, 'mitoinfo')
    end

    %clearvars -except mitoWidth mitoArea mitoLength mitoAR mitoOXPHOS mitoOMP mitoOXPHOSparam
end
