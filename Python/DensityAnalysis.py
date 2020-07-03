import os
import cv2 as cv
import glob
import tifffile
import read_roi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tkinter.filedialog import askdirectory
from skimage import measure

dirpath = askdirectory(title='Choose your folder...',initialdir='X:/Mitography/NEW/Antimycin Treatments_April2020/6h_5nM AA')  # directory path
files_neur = glob.glob(os.path.join(dirpath,'Image_0*NeuritesBinary.tif'))
files_neurlen = glob.glob(os.path.join(dirpath,'Image_0*NeuritesLength.txt'))
files_mitoana = glob.glob(os.path.join(dirpath,'Image_0*MitoAnalysisFull.txt'))
files_oxphos = glob.glob(os.path.join(dirpath,'*OXPHOSbool.xlsx'))
files_mitobin = glob.glob(os.path.join(dirpath,'*MitoBinary.tif'))
files_roi = glob.glob(os.path.join(dirpath,'*RoiSet.zip'))

for filepath_neur in files_neur:
    # create empty dataframe for all the mito of interest to be saved in
    mitoofinterest = pd.DataFrame()

    # print imgname, get index of img in all lists, and start loading data
    imgname_neur = filepath_neur.split('\\')[1].split('.')[0]
    print(imgname_neur)
    imgno = int(imgname_neur.split('_')[1].split('-')[0])
    imgidx = files_neur.index(filepath_neur)

    # read tot neurite length and pixel size
    [throw, pxs_nm, totneurlen] = np.loadtxt(files_neurlen[imgidx])

    # read binary mitochondria image and get labelled binary mito image
    with tifffile.TiffFile(files_mitobin[imgidx]) as tif:
        mitobin = tif.pages[0].asarray().astype('uint8')  # image as numpy array
        labelmito = measure.label(mitobin)

    # read ROIs of mitochondria of interest
    rois = read_roi.read_roi_zip(files_roi[imgidx])
    # get the coordinates from ~the middle of the ROIs
    xcoor_rois = []; ycoor_rois = []; mitonums = []
    for roi in rois:
        roi_xs = rois[roi]['x']
        roi_ys = rois[roi]['y']
        xroi = (roi_xs[0]+roi_xs[len(roi_xs)//2])//2
        yroi = (roi_ys[0]+roi_ys[len(roi_ys)//2])//2
        xcoor_rois.append(xroi)
        ycoor_rois.append(yroi)
    # compare coords with labelled mito img (labelled in same order as MitoBinary), and save mito numbers to list
    for x,y in zip(xcoor_rois,ycoor_rois):
        mitonum = xcoor_rois.index(x)
        mitonums.append(labelmito[y,x])
    # save mitonums of interest to all-data dataframe
    mitoofinterest['mitonum'] = mitonums

    # read OXPHOS boolean data
    oxphos = pd.read_excel(files_oxphos[imgidx], sheet_name='OXPHOS')
    oxphos = list(oxphos.iloc[:,-1][:len(rois)])
    # save OXPHOS data to all-data dataframe
    mitoofinterest['OXPHOSbool'] = oxphos

    # read the mitoanalysis file and get coordinates and area for all mitos in img
    anafull = pd.read_csv(files_mitoana[imgidx], sep='\t', header=None)
    xcor_px = round(anafull[0]/pxs_nm*1000).astype(int)
    ycor_px = round(anafull[1]/pxs_nm*1000).astype(int)
    area = anafull[3]
    # gather the area of the mito of interest in a list
    mitoareas = []
    for mitono in mitoofinterest['mitonum']:
        mitoareas.append(area[mitono-1])
    # add list of areas to all-data dataframe
    mitoofinterest['area'] = mitoareas

    # read binary neurites image
    with tifffile.TiffFile(files_neur[imgidx]) as tif:
        neurbin = tif.pages[0].asarray().astype('uint8')  # image as numpy array
    # mark in a list which mito are in mitoneur and which are not, of the mitos of interest
    mitoinneur = []
    for mitono in mitoofinterest['mitonum']:
        mitoneur = neurbin[ycor_px[mitono-1]][xcor_px[mitono-1]]
        mitoinneur.append(mitoneur)
    # add list of mitoneur to all-data dataframe
    mitoofinterest['inneurite'] = mitoinneur

    # save all-data dataframe to csv
    mitoofinterest.to_csv(os.path.join(dirpath,f'Image_{imgno:03d}-MitoOfInterest.csv'), header=True, index=False)


    ## save which mito are in neurites
    #pd.DataFrame(mitoinneur).to_csv(os.path.join(dirpath,f'Image_{imgno:03d}-MitoInNeurites.txt'), header=False, index=False)

    """
    # calculate number of mitos and small mitos per neurite length
    areathresh = 0.086
    smallmito = np.array(area)<areathresh
    mitoperum = np.sum(mitoinneur)/totneurlen
    smallmitoperum = np.sum(smallmito & list(map(bool,mitoinneur)))/totneurlen

    # save results
    np.savetxt(os.path.join(dirpath,f'Image_{imgno:03d}-Density.txt'), [mitoperum, areathresh, smallmitoperum], fmt='%.4f')
    """

    """
    # visualize and mark in a list which small mito are in mitoneur and which are not
    marker = ['ro','go','bo']
    fig = plt.figure(figsize=(10,10))
    plt.imshow(neurbin)
    for mitono in range(0,len(xcor_px)):
        if mitoinneur[mitono]:
            if smallmito[mitono]:
                plt.plot(xcor_px[mitono], ycor_px[mitono], marker[1], markersize=10)
            else:
                plt.plot(xcor_px[mitono], ycor_px[mitono], marker[2], markersize=10)
        else:
            plt.plot(xcor_px[mitono], ycor_px[mitono], marker[0], markersize=10)
    """
#mitoperum, smallmitoperum

