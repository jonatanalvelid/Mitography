%%%
% Mitography - TFAM analysis
% Analysing number of nucleoids in mitochondria, based on TFAM and
% OMM-staining.
%----------------------------
% Version: 191203
% Last updated features: New script
%
% @jonatanalvelid
%%%

clear

% Add function folder to filepath, so that those functions can be read.
functionFolder = fileparts(which('findFunctionFolders.m'));
addpath(genpath(functionFolder));

masterFolderPath = strcat(uigetdir('X:\Mitography'),'\');
fileList = dir(fullfile(masterFolderPath, 'Image_*.txt'));
for i = 1:length(fileList)
    filenumbers(i) = str2num(fileList(i).name(7:9));
end
lastFileNumber = max(filenumbers);

threshsize = 8;  % Lower threshold size in pixels for binary mitochondria

filenameallPxs = '_PixelSizes.txt';
filenameallMito = '_MitoAnalysis.txt';
filenamenucleoids = '_Nucleoids.txt';
filenameAnalysisSave = '_MitoAnalysisFull.txt';
filenameMitoBinary = '_MitoBinary.tif';
filenameSomaBinary = '-SomaBinary.tif';
filenameAISBinary = '_AISBinary.tif';
filenameAxonDistInfo = '_AxonDistInfo.txt';
fileNumbers = 1:lastFileNumber;

for fileNum = fileNumbers
    filepathpxs = strFilepath(fileNum,filenameallPxs,masterFolderPath);
    filepathmito = strFilepath(fileNum,filenameallMito,masterFolderPath);
    filepathnucleoids = strFilepath(fileNum,filenamenucleoids,masterFolderPath);
    filepathMitoBinary = strFilepath(fileNum,filenameMitoBinary,masterFolderPath);
    filepathAnaSave = strFilepath(fileNum,filenameAnalysisSave,masterFolderPath);
    filepathSomaBinary = strFilepath(fileNum,filenameSomaBinary,masterFolderPath);
    filepathAISBinary = strFilepath(fileNum,filenameAISBinary,masterFolderPath);
    filepathAxonDistInfo = strFilepath(fileNum,filenameAxonDistInfo,masterFolderPath);
    
    try
        % Read the mito and line profile data
        datamito = dlmread(filepathmito,'',1,1);
        datanucleoids = dlmread(filepathnucleoids,'',1,1);
        [~,params] = size(datamito);
        
        % Read the pixel size (in nm)
        datapxs = dlmread(filepathpxs,'',1,1);
        pixelsize = datapxs(1,1)/1000;
        
        % Load binary mitochondria image 
        imagemitobinary = imread(filepathMitoBinary);
        imsize = size(imagemitobinary);
        
        % Make binary nucleoid center map 
        nucleoidmap = zeros(size(imagemitobinary));
        [numnucl, ~] = size(datanucleoids);
        for i = 1:numnucl
            xpos = datanucleoids(i,2);
            xpos = round(xpos/pixelsize);
            ypos = datanucleoids(i,3);
            ypos = round(ypos/pixelsize);
            % Make sure all coordinates are in the range of the img size
            xpos = min(max(xpos,1),imsize(2));
            ypos = min(max(ypos,1),imsize(1));
            nucleoidmap(ypos,xpos) = 1;
        end
        
        % Remove small objects and make labelled binary mitochondria image
        imagemitobinary = bwareaopen(imagemitobinary, threshsize);
        [labelmito, num] = bwlabel(imagemitobinary');
        labelmito = labelmito';
        
        % Mark those mitochondria that are in "soma" areas
        try
            imagesomabinary = imread(filepathSomaBinary);
            imagesomabinary = logical(imagesomabinary);
        catch err
            imagesomabinary = zeros(size(imagemitobinary));
            imagesomabinary = logical(imagesomabinary);
        end
        for i = 1:num
            xpos = datamito(i,1);
            xpos = round(xpos/pixelsize);
            ypos = datamito(i,2);
            ypos = round(ypos/pixelsize);
            % Make sure all coordinates are in the range of the img size
            xpos = min(max(xpos,1),imsize(2));
            ypos = min(max(ypos,1),imsize(1));
            datamito(i,params+1) = imagesomabinary(ypos,xpos);
        end

        % Get the distance from the soma along the axon to the mitochondria
        %%% MANUALLY CREATE .txt file that carries information about the
        %%% seed point for the bwdistgeodesic transformation, and the
        %%% previous distance along the axon.
        % Read the axon distance info
        datadistinfo = dlmread(filepathAxonDistInfo,'',1,1);
        seedx = datadistinfo(1); seedy = datadistinfo(2); prevdist = datadistinfo(3);
        % Read the binary AIS-image (axon image)
        imageaisbin = imread(filepathAISBinary);
        imageaisbin = logical(imageaisbin);
        aisdist = bwdistgeodesic(imageaisbin, seedx, seedy, 'quasi-euclidean') + prevdist;
        for i = 1:num
            xpos = datamito(i,1);
            xpos = round(xpos/pixelsize);
            ypos = datamito(i,2);
            ypos = round(ypos/pixelsize);
            % Make sure all coordinates are in the range of the img size
            xpos = min(max(xpos,1),imsize(2));
            ypos = min(max(ypos,1),imsize(1));
            datamito(i,params+2) = aisdist(ypos,xpos) * pixelsize;
        end
        % Round the distances to three decimals
        datamito(:,8) = round(datamito(:,8),3);
        
        % Get number of nucleoids in each mito and save to mitoinfo
        for i = 1:num
            singlemitobinary = ismember(labelmito, i);
            singlemitonucleoids = nucleoidmap.*singlemitobinary;
            datamito(i,params+3) = sum(sum(singlemitonucleoids));
        end

        % Save data
        disp(strcat(num2str(fileNum),': Done.'))
%         dlmwrite(filepathAnaSave,datamito,'delimiter','\t');
        writematrix(datamito,filepathAnaSave,'Delimiter','tab')
    catch err
        disp(strcat(num2str(fileNum),': General error.'));
    end 
    
end