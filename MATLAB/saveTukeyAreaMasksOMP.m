%%%
%
% Mitography - PEX (AA/control) - Tukey bkg calc mask for mitomarker
%
%----------------------------
%
% Version: 210301
%
%%%

clear

% Add functions folder to filepath and get data folder path
filename = matlab.desktop.editor.getActiveFilename;
parentfolder = getfield(fliplr(regexp(fileparts(fileparts(filename)),'/','split')),{1});
doubleparentfolder = getfield(fliplr(regexp(fileparts(fileparts(fileparts(filename))),'/','split')),{1});
functionsfolder = fullfile(parentfolder{1},'functions');
addpath(functionsfolder);

masterFolderPath = strcat(uigetdir('D:\Data analysis\Mitography\AA-PEX\analysis-pex_start\exp2ct\pex-matlab_analysis'),'\');

lastFileNumber = input('What is the number of the last image? ');
mitosPerFile = 1000;

imgNumbers = 1:lastFileNumber;

%% 
%%% ANALYSIS

filenamePEXBinary = '-pexbinary.tif';
filenameNeuritesBinary = '-neuritesbinary.tif';
filenameMitoBinary = '-mitobinary.tif';
filenameSave = '-mitobkgbinary.tif';

fileNum = 1;
for imgNum = imgNumbers
    filepathPexBinary = strFilepath(imgNum,filenamePEXBinary,masterFolderPath);
    filepathNeuritesBinary = strFilepath(imgNum,filenameNeuritesBinary,masterFolderPath);
    filepathMitoBinary = strFilepath(imgNum,filenameMitoBinary,masterFolderPath);
    savepath = strFilepath(imgNum,filenameSave,masterFolderPath);
 
    try
        % Read images
        imagepexbinary = imread(filepathPexBinary);
        imagepexbinary = logical(imagepexbinary);
        try
            imageneuritesbinary = imread(filepathNeuritesBinary);
            imageneuritesbinary = logical(imageneuritesbinary);
        catch err
            imageneuritesbinary = zeros(size(imagepexbinary));
            imageneuritesbinary = logical(imageneuritesbinary);
        end
        try
            imagemitobinary = imread(filepathMitoBinary);
            imagemitobinary = logical(imagemitobinary);
        catch err
            imagemitobinary = zeros(size(imagepexbinary));
            imagemitobinary = logical(imagemitobinary);
        end

        bkginsideMask = imageneuritesbinary & ~imagepexbinary & ~imagemitobinary;
        
        imwrite(bkginsideMask, savepath)
        
    catch err
        disp(strcat(num2str(imgNum),': No image with this number or a file reading error.'))
    end
end

