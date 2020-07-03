import os
import cv2 as cv
import glob
import tifffile
import numpy as np
import pandas as pd
from tkinter.filedialog import askdirectory
import matplotlib.pyplot as plt

dirpath = askdirectory(title='Choose your folder...',initialdir='X:/Mitography/NEW/Antimycin Treatments_April2020/6h_5nM AA')  # directory path
files_neur = glob.glob(os.path.join(dirpath,'Image_0*NeuritesBinary.tif'))
files_neurlen = glob.glob(os.path.join(dirpath,'Image_0*NeuritesLength.txt'))
files_mitoana = glob.glob(os.path.join(dirpath,'Image_0*MitoAnalysisFull.txt'))

for filepath_neur in files_neur:
    imgname_neur = filepath_neur.split('\\')[1].split('.')[0]
    print(imgname_neur)
    imgno = int(imgname_neur.split('_')[1].split('-')[0])
    imgidx = files_neur.index(filepath_neur)
    # read tot neurite length and pixel size
    [throw, pxs_nm, totneurlen] = np.loadtxt(files_neurlen[imgidx])

    # read binary neurites image
    with tifffile.TiffFile(files_neur[imgidx]) as tif:
        neurbin = tif.pages[0].asarray().astype('uint8')  # image as numpy array

    # read the mitoanalysis file and get coordinates and area for mitos
    anafull = pd.read_csv(files_mitoana[imgidx], sep='\t', header=None)
    xcor_px = round(anafull[0]/pxs_nm*1000).astype(int); ycor_px = round(anafull[1]/pxs_nm*1000).astype(int); area = anafull[3]

    # visualize and mark in a list which mito are in mitoneur and which are not
    mitoinneur = []
    for mitono in range(0,len(xcor_px)):
        mitoneur = neurbin[ycor_px[mitono]][xcor_px[mitono]]
        mitoinneur.append(mitoneur)

    # save which mito are in neurites
    pd.DataFrame(mitoinneur).to_csv(os.path.join(dirpath,f'Image_{imgno:03d}-MitoInNeurites.txt'), header=False, index=False)

    # calculate number of mitos and small mitos per neurite length
    areathresh = 0.086
    smallmito = np.array(area)<areathresh
    mitoperum = np.sum(mitoinneur)/totneurlen
    smallmitoperum = np.sum(smallmito & list(map(bool,mitoinneur)))/totneurlen

    # save results
    np.savetxt(os.path.join(dirpath,f'Image_{imgno:03d}-Density.txt'), [mitoperum, areathresh, smallmitoperum], fmt='%.4f')


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

