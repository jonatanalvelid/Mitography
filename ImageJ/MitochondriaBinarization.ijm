// Binarize STED mitochondria image

noImages = 1;
outerproffactor = 3;
//seems to be necessary in fixed samples where the membrane labeling is more dotty (OMP25 vs Tom20 rather than live vs fixed).
omp25 = 1;  //only one of these two should be 1, the other 0
tom20 = 0;  //only one of these two should be 1, the other 0
allmito = 1;  //for counting all mito, or deleting the smallest and big ones for example, to get less noise
localthreshold = 10;  //try local thresholding, Bernsen variant. Seems to work nicely on well-labelled OMP25 images.
localthresholdmethod = "Bernsen";  //local thresholding method to use.

getPixelSize(unit, pixelWidth, pixelHeight);
getDimensions(width, height, channels, slices, frames);

run("Set Scale...", "distance=1 known="+pixelWidth+" pixel=1 unit=micron");
	
wait(500);

imnameor = getTitle();
imname = substring(imnameor, 0, 9);

//Taking the images and renaming them in a standard way.
rename("MitoOriginalImageSoma");

//Starting the analysis by binarizing the mitochondria image and calculating important parameters for them.
selectWindow("MitoOriginalImageSoma");
run("Duplicate...", "title=MitoOriginalImageSoma");
rename("mitobinaryalt");
selectWindow("mitobinaryalt");
if(tom20 == 1) {
	run("Gaussian Blur...", "sigma=0.05 scaled");
	if(localthreshold == 1) {
		run("Auto Local Threshold", "method=" + localthresholdmethod + " radius=15 parameter_1=0 parameter_2=0 white");
	} else {
		run("Make Binary");
		run("Fill Holes");
	}
} else if(omp25 == 1) {
	run("Gaussian Blur...", "sigma=0.1 scaled"); //Try this, for the OMP25 labeling at least
	if(localthreshold == 1) {
		run("Auto Local Threshold", "method=" + localthresholdmethod + " radius=15 parameter_1=0 parameter_2=0 white");
	} else {
		run("Make Binary");
		run("Fill Holes");
	}
}
selectWindow("mitobinaryalt");
rename("mitobinary");
selectWindow("mitobinary");

if(allmito == 1) {
	if(tom20 == 1) {
		run("Analyze Particles...", "size=0.015-Infinity display clear include add"); //for counting all, big and small MDVs
	} else if(omp25 == 1) {
		run("Analyze Particles...", "size=0.003-Infinity display clear include add"); //for counting all, big and small MDVs
	}
} else if(allmito == 0) {
	//run("Analyze Particles...", "size=0.08-3.00 display exclude clear include add"); - for normal
	//run("Analyze Particles...", "size=0.08-Infinity display exclude clear include add"); //for getting all the mitochondria areas, also the big ones
	//run("Analyze Particles...", "size=0.02-Infinity display exclude clear include add"); //for getting all the mitochondria areas, also the big and small ones
	run("Analyze Particles...", "size=0.08-Infinity display exclude clear include add"); //for confocal
}
updateResults; 

selectWindow("mitobinary");
run("Duplicate...", "title=mitobinary2");
selectWindow("mitobinary2");
if(allmito == 1) {
	if(tom20 == 1) {
		run("Analyze Particles...", "size=0.015-Infinity display clear include add"); //for counting all, big and small MDVs
	} else if(omp25 == 1) {
		run("Analyze Particles...", "size=0.003-Infinity display clear include add"); //for counting all, big and small MDVs
	}
} else if(allmito == 0) {
	//run("Analyze Particles...", "size=0.08-3.00 display exclude clear include add"); - for normal
	//run("Analyze Particles...", "size=0.08-Infinity display exclude clear include add"); //for getting all the mitochondria areas, also the big ones
	//run("Analyze Particles...", "size=0.02-Infinity display exclude clear include add"); //for getting all the mitochondria areas, also the big and small ones
	run("Analyze Particles...", "size=0.08-Infinity display exclude clear include add"); //for confocal
}
run("Flatten");

length1 = nResults();

wait(200);
selectWindow("mitobinary");

run("Close");
run("Close");
