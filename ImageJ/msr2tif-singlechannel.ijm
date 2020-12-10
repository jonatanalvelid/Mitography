dir = getDirectory("Choose the directory");
savedir = "E:\\PhD\\data_analysis\\Temp\\";
filelist = getFileList(dir);
//Array.sort(filelist);
filenamebase = "\\"+dir+"\\";

function leftPad(s, width) {
  while (lengthOf(s)<width) {
      s = "0"+s;
  }
  return toString(s);
}

function saveimages(filename, savedir, savenum, namesuffix) {
		run("Bio-Formats Importer", "open=" + filename + " color_mode=Default view=Hyperstack stack_order=XYCZT series_2");
		savenumtext = leftPad(savenum,3);
		savename="Image_"+savenumtext+"-"+namesuffix+".tif";
		saveAs("Tiff", savedir+savename);
		close();
}

for(r=0;r<filelist.length;r++) {
	if(endsWith(filelist[r],'.msr')) {
		print(filelist[r]);
		filepath = filenamebase+filelist[r];
		print(filepath);
		filename_split = split(filelist[r],'-');
		filesavenum = filename_split[1];
		filesavenum = split(filesavenum,'.');
		filesavenum = filesavenum[0];
		filenameprefix = filename_split[0];
		filenameprefix = toString(filenameprefix);
		saveimages(filepath, savedir, ""+filesavenum, "mitochondria");
	}
}
