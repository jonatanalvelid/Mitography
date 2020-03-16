%%%
%
% Mitography main code. Analyses all info and line profiles spitted out
% from the ImageJ Mitography script, and gives the info per mitochonrdia
% out.
%
%----------------------------
%
% Version: 200313
% New script to handle TMR data, with TMR binary maps, with a boolean
% output variable saying if the mitochondria has TMR or not, and the
% average TMR signal inside the mitochondria.
%
%%%

clear

% Add function folder to filepath, so that those functions can be read.
functionFolder = fileparts(which('findFunctionFolders.m'));
addpath(genpath(functionFolder));

masterFolderPath = strcat(uigetdir('X:\Mitography\NEW\Metabolism\Control\OMP25'),'\');
filenameparam = 'ImageJAnalysisParameters.txt';
filepathparam = strcat(masterFolderPath,filenameparam);
try
    dataparam = readtable(filepathparam);
    mitoLineProfLength = dataparam.MitoLineProfLen;
    actinLineProfLength = dataparam.ActLineProfLen;
catch err
    mitoLineProfLength = input('What is the mitochondria line profile length (um)? ');
    actinLineProfLength = input('What is the actin line profile length (um)? ');
end

% lastFileNumber = input('What is the number of the last image? ');
mitosPerFile = 1000;

filenameallPxs = '_PixelSizes.txt';
filenameAnalysis = '_MitoAnalysis.txt';
filenameAnalysisSave = '_MitoAnalysisFull.txt';
filenameMitoBinary = '_MitoBinary.tif';
filenameSomaBinary = '-SomaBinary.tif';
filenameTMR = '-TMR.tif';

for fileNum = 1%fileNumbers
    
    filepathAnaSave = strFilepath(fileNum,filenameAnalysisSave,masterFolderPath);
    filepathAna = strFilepath(fileNum,filenameAnalysis,masterFolderPath);
    filepathpxs = strFilepath(fileNum,filenameallPxs,masterFolderPath);
    filepathTMR = strFilepath(fileNum,filenameTMR,masterFolderPath);
    filepathMitoBinary = strFilepath(fileNum,filenameMitoBinary,masterFolderPath);
    filepathSomaBinary = strFilepath(fileNum,filenameSomaBinary,masterFolderPath);
    
    try
        try
            dataAnalysis = dlmread(filepathAna,'',1,1);
            % dataMitoNoStr = dlmread(filepathMitoNonStr,'',1,1);
            sizeData = size(dataAnalysis);
        catch err
            disp(strcat(num2str(fileNum),': File reading error.'));
        end
        % Read the pixel size
        datapxs = dlmread(filepathpxs,'',1,1);
        pixelsize = datapxs(1,1);
%         pixelsize = 0.0292969;

        
        %%% MITOCHONDRIA TMR CHECK AND FLAGGING
        
        % Load TMR image
        imagetmr = imread(filepathTMR);
        % Load binary mitochondria image 
        imagemitobinary = imread(filepathMitoBinary);
        imsize = size(imagemitobinary);
        % Remove small objects and make labelled binary mitochondria image
        % flip axes of imagemitobinary so numbering matches previous
        % numbering in imageJ, i.e. in mitoAnalysis.text file
        [labelmito, num] = bwlabel(imagemitobinary');
        % flip axes back
        labelmito = labelmito';
        % Get single mito map, multiply with tmr image
        for i = 1:sizeData(1)
        end

        if not(isempty(dataAnalysis))
            for i=1:sizeData(1)
                % get binary img of single mitochondria
                singlemitobinary = ismember(labelmito, i);
                % get a list of tmr pixels in this area
                tmrsignal = imagetmr(singlemitobinary);
                % get average TMR signal/pixel per mito, save to data (110)
                tmrsignalavg = mean(tmrsignal);
                dataAnalysis(i,10) = tmrsignalavg;
%                 xcor = max(min(round(dataAnalysis(i,1)/pixelsize),imsize(1)),1);
%                 ycor = max(min(round(dataAnalysis(i,2)/pixelsize),imsize(1)),1);
%                 dataAnalysis(i,12) = singlemitobinary(ycor,xcor);
            end
            
            %decide what is the limiting TMR signal for positive/negative
            %signal, and save another variable (111) - boolean yes/no
            threshsignal = 0.5;  % decide thresh signal here
            
            % find thresh by double gaussian fit and find intersection %%%
            % FIX THIS GAUSSIAN FIT HERE NEXT TIME
            numgroups = 15;
            [cnts,edges] = histcounts(dataAnalysis(:,10),numgroups);
            edges(end) = [];
            cnts = cnts/max(cnts);
            x = edges + (edges(2)-edges(1))/2;
            [xData, yData] = prepareCurveData( x, cnts );

            % Set up fittype and options.
            ft = fittype( 'a1*exp(-((-x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)+d', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 0 -x(end)/4 0 -5 -5 0];
            opts.StartPoint = [1 0.5 0 x(round(numgroups/2)) -0.5 3 0.1];
            opts.Upper = [10 5 x(end)/4 x(end) 10 10 0.3];

            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            cfs = coeffvalues(fitresult);
            disp(cfs)
            
            xsampl = 1:x(end)/100:x(end);
            bkggauss = cfs(1).*exp(-((-xsampl-cfs(3))./cfs(5)).^2);
            signalgauss = cfs(2).*exp(-((-xsampl-cfs(4))./cfs(6)).^2);
            [~,idx] = min(bkggauss-signalgauss);
            disp(xsampl(idx))
            
            % Plot fit with data.
            figure( 'Name', 'untitled fit 1' );
            h = plot( fitresult, xData, yData );
            legend( h, 'cnts vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
            % Label axes
            xlabel( 'x', 'Interpreter', 'none' );
            ylabel( 'cnts', 'Interpreter', 'none' );
            grid on
           
            
            for i=1:sizeData(1)
                if dataAnalysis(i,10) > threshsignal
                    dataAnalysis(i,11) = 1;
                else
                    dataAnalysis(i,11) = 0;
                end
            end
        end
        
%         disp(strcat(num2str(fileNum),': Data handling done.'))
%         dlmwrite(filepathAnaSave,dataAnalysis,'delimiter','\t');
    catch err
        disp(strcat(num2str(fileNum),': No image with this number or a file reading error.'))
    end
end