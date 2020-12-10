clear all

dir = strcat('E:\PhD\data_analysis\temp_pex\raw_results','\');
file = 'pex_e2_ct_c4.mat';
load(strcat(dir,file))

mitoinfo = struct;
mitoinfo.exp = 2;
mitoinfo.ctaa = 'ct';
mitoinfo.cell = 4;
mitoinfo.ar = mitoAR;
mitoinfo.area = mitoArea;
mitoinfo.length = mitoLength;
mitoinfo.omp = mitoOMP;
mitoinfo.pex = mitoPEX;
mitoinfo.pexparam = mitoPEXparam;
mitoinfo.width = mitoWidth;

save(file, 'mitoinfo')