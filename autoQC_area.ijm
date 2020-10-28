// autoQC_area.ijm
// - enhancement to autoQC_counter.ijm for single images only, displaying area masks 
//
// Usage:
// - single image and batch modes, selected in dialog
// - user-specified red-channel (mean +/- N stdDev, for peaks) & ratio thresholds 
// - both modes expect 2D multi-channel (2+) input image(s)
// - ROIs for cells *only* should be in ROI manager, or image overlays in batch mode
// - batch mode runs on a folder of .tif images saved from ImageJ with ROI overlays
//
// Output:
// - output image with cell ROIs and multi-point spot ROIs added to overlay:
//   - with additional "redOnlySpotMask" and "redAndYellowSpotMask" output channels
//   - output channel 1=red, 2=green, 3=red+yellow (yellow), 4=red only (magenta)
// - summary results: 1 row per cell ROI (results for all images appended to table):
//   - Image name
//   - Cell ROI name
//   - ROI XY top-left coords, X_tl & Y_tl
//   - cell area
//   - number of spots below red threshold intensity, nBlack
//   - number of red spots, nRed
//   - number of red+green spots, nYellow
//   - total areas of: red, yellow, red+yellow
//   - (parameters: peakNstd,redMinimumSD,rThresh,spotSize)
//
// g.ball@dundee.ac.uk, Dundee Imaging Facility (2020)
// license: MIT License



// parameters
grnCh = 1;  // channel for grnChondrial marker (green)
redCh = 2;  // channel for mitophagy marker (red)
peakNstd = 2;  // standard deviations above mean for peaks (red channel)
redMinimumSD = 2.5;  // minimum red intensity for spots: mean+/-nStdDev
rThresh = 1.0;  // minimum red/green value for red spots
spotSize = 0.1;  // calibrated units (microns)
DoGmultiplier = 5;  // multiplier for DoG filter
showRedSpotImage = true;  // option to show red spot image


// --- start Macro ---

// check we have at least 1 ROI & composite image before going further
nRois = roiManager("count");
if (nRois < 1) {
	exit("Please add at least 1 cell outline to ROI Manager!");
}
if (!is("composite")) {
	exit("Macro requires multi-channel image stack!");
}

run("Select None");
imageTitle = getTitle();
Stack.getDimensions(width, height, nc, nz, nt);
if (nc < 2) {
	exit("Macro requires at least 2 channels! (red & green)");
}
inpID = getImageID();  // input image
inpTitle = getTitle();
Stack.getUnits(xu, yu, zu, tu, vu);
areaUnits = "" + xu + "2";
inpBitDepth = bitDepth();


// dialog to check/update channels & options
Dialog.create("autoQC_area");
Dialog.addNumber("Green channel", grnCh);
Dialog.addNumber("Red channel", redCh);
Dialog.addNumber("Spot size in calibrated units (microns?)", spotSize);
Dialog.addNumber("Number of stdDevs above background for peaks", peakNstd);
Dialog.addNumber("Minimum red intensity for spots (stdDevs above mean) ", redMinimumSD);
Dialog.addNumber("Minimum red/green intensity ratio for red spots", rThresh);
Dialog.show();
grnCh = Dialog.getNumber();
redCh = Dialog.getNumber();
spotSize = Dialog.getNumber();
peakNstd = Dialog.getNumber();
redMinimumSD = Dialog.getNumber();
rThresh = Dialog.getNumber();


setBatchMode("hide");  // hide images during calculations
setOption("BlackBackground", true);

// pre-process with median filter
run("Select None");
run("Duplicate...", "duplicate");
filteredID = getImageID();
run("Median...", "radius=1");  // all channels
filteredTitle = "" + baseName(inpTitle) + "_filtered.tif";
rename(filteredTitle);

// generate red spot image (DoG filter)
Stack.setChannel(redCh);
redSpotImageID = dogFilter(filteredID, spotSize, DoGmultiplier);
rename("redSpotImage");

run("Clear Results");
resultsRow = 0;
// iterate over ROIs (cells)
for (r = 0; r < nRois; r++) {
	// find peaks in DoG-filtered red "spot image"
	selectImage(redSpotImageID);
	run("Select None");
	getStatistics(area, mean, min, max, std);  // whole image
	peakThresh = peakNstd * std;  
	roiManager("select", r);
	run("Find Maxima...", "prominence=" + peakThresh + " output=[Point Selection]");
	getSelectionCoordinates(xc, yc);
	nSpots = xc.length;
	// create red spot image mask
	run("Select None");
	run("Duplicate...", " ");
	setThreshold(peakThresh, max);
	run("Convert to Mask");
	rename("redSpotMask");
	redSpotMaskID = getImageID();
	// create red/green ratio image mask
	selectImage(filteredID);
	Stack.setChannel(redCh);
	run("Duplicate...", " ");
	rename("redChannelImage");
	redImageID = getImageID();
	selectImage(filteredID);
	Stack.setChannel(grnCh);
	run("Duplicate...", " ");
	rename("greenChannelImage");
	greenImageID = getImageID();
	imageCalculator("Divide create 32-bit", "redChannelImage","greenChannelImage");
	getStatistics(area, rMean, rMin, rMax, rStd);
	setThreshold(rThresh, rMax);
	run("Convert to Mask");
	rename("redOverGreenRatioMask");
	ratioMaskID = getImageID();
	selectImage(greenImageID);
	close();
	// calculate redMinimum using whole-image mean, SD; & create redMinimumMask
	selectImage(redImageID);
	run("Select None");
	getStatistics(area, redMean, redMin, redMax, redStd);  // whole image
	redMinimum = redMean + (redStd * redMinimumSD);
	setThreshold(redMinimum, max);
	run("Convert to Mask");
	rename("redAboveMinimumMask");
	redAboveMinimumID = getImageID();
	// calculate mask for above-threshold red spots (i.e. red & yellow spots)
	imageCalculator("AND create", "redSpotMask","redAboveMinimumMask");
	rename("redAndYellowSpotMask");
	redAndYellowSpotMaskID = getImageID();
	selectImage(redSpotMaskID);
	close();
	selectImage(redAboveMinimumID);
	close();
	areaRedAndYellow = measureMaskArea(redAndYellowSpotMaskID, r);
	selectImage(redAndYellowSpotMaskID);
	run("" + inpBitDepth + "-bit");  // convert to input bit-depth for merge
	// calculate mask for above-thresh red spots that are red (R/G above ratio thresh)
	imageCalculator("AND create", "redAndYellowSpotMask","redOverGreenRatioMask");
	rename("redOnlySpotMask");
	redOnlyMaskID = getImageID();
	selectImage(ratioMaskID);
	close();
	areaRedOnly = measureMaskArea(redOnlyMaskID, r);
	selectImage(redOnlyMaskID);
	run("" + inpBitDepth + "-bit");  // convert to input bit-depth for merge
	// append red+yellow and red-only mask channels to filtered image
	selectImage(filteredID);
	run("Split Channels");
	C1name = "c1=[C" + redCh + "-" + filteredTitle + "]";
	C2name = "c2=[C" + grnCh + "-" + filteredTitle + "]";
	C3name = "c3=redAndYellowSpotMask";
	C4name = "c4=redOnlySpotMask";
	channelsString = "" + C1name + " " + C2name + " " +  C3name + " " + C4name + " create";
	run("Merge Channels...", channelsString);
	filteredID = getImageID();
	rename(filteredTitle);
	Stack.setChannel(1);
	run("Red");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(2);
	run("Green");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(3);
	run("Yellow");
	setMinAndMax(0, 255);
	Stack.setChannel(4);
	run("Magenta");
	setMinAndMax(0, 255);
	// record cell ROI stats
	selectImage(filteredID);
	roiManager("select", r);
	Overlay.addSelection();
	setResult("CellROI", resultsRow, Roi.getName());
	Roi.getBounds(x, y, w, h);
	setResult("X_tl", resultsRow, "" + x);
	setResult("Y_tl", resultsRow, "" + y);
	getStatistics(area, mean, min, max, std);  // cell ROI
	setResult("CellArea_" + areaUnits, resultsRow, area);
	nBlack = 0;
	xBlack = newArray(nSpots);
	yBlack = newArray(nSpots);
	nRed = 0;
	xRed = newArray(nSpots);
	yRed = newArray(nSpots);
	nYellow = 0;
	xYellow = newArray(nSpots);
	yYellow = newArray(nSpots);
	// get green & red intensities (Stack.setChannel() is slow...)
	iGrn = newArray(nSpots);
	Stack.setChannel(2);  // green (post-merge)
	for (i = 0; i < xc.length; i++) {
		iGrn[i] = getPixel(xc[i], yc[i]);
	}
	iRed = newArray(nSpots);
	Stack.setChannel(1);  // red (post-merge)
	for (i = 0; i < xc.length; i++) {
		iRed[i] = getPixel(xc[i], yc[i]);
	}
	// iterate over points classifying
	for (i = 0; i < xc.length; i++) {
		if (iRed[i] < redMinimum) {
			xBlack[nBlack] = xc[i];
			yBlack[nBlack] = yc[i];
			nBlack++;
		} else {
			if (iRed[i] > iGrn[i]*rThresh) {
				xRed[nRed] = xc[i];
				yRed[nRed] = yc[i];
				nRed++;
			} else {
				xYellow[nYellow] = xc[i];
				yYellow[nYellow] = yc[i];
				nYellow++;
			}
		}
	}
	// write stats row to results table
	setResult("nBlack", resultsRow, nBlack);
	setResult("nRed", resultsRow, nRed);
	setResult("nYellow", resultsRow, nYellow);
	setResult("areaRed_" + areaUnits, resultsRow, areaRedOnly);
	setResult("areaRedAndYellow_" + areaUnits, resultsRow, areaRedAndYellow);
	setResult("areaYellow_" + areaUnits, resultsRow, areaRedAndYellow - areaRedOnly);	
	setResult("spotSize_" + areaUnits, resultsRow, spotSize);
	setResult("peakNstd", resultsRow, peakNstd);
	setResult("rThresh", resultsRow, rThresh);
	setResult("redMinimumSD", resultsRow, redMinimumSD);
	setResult("redMinimum", resultsRow, redMinimum);
	updateResults();
	resultsRow++;
	// Add coloured classified spots to overlay
	addSpotsToOverlay(xBlack, yBlack, nBlack, "black");
	addSpotsToOverlay(xRed, yRed, nRed, "red");
	addSpotsToOverlay(xYellow, yYellow, nYellow, "yellow");
	run("Select None");
}
selectImage(redSpotImageID);
if (showRedSpotImage) {
	run("Select None");
	run("Enhance Contrast", "saturated=0.35");
} else {
	close();
}
paramString = "R" + redCh + "_G" + grnCh;
IJ.renameResults("Results-autoQC-" + paramString);

setBatchMode("exit and display");


// --- function definitions ---

function dogFilter(inputImageID, scale, multiplier) {
	// return imageID of new, DoG-filtered image for current channel
	// - "scale" of smallest gaussian in calibrated units
	// - N.B. does "Select None" before processing
	selectImage(inputImageID);
	outputTitle = "DoG" + scale + "x" + multiplier + "_" + getTitle();
	run("Select None");
	run("Duplicate...", " ");
	run("32-bit");
	baseImageID = getImageID;
	run("Duplicate...", " ");
	run("Gaussian Blur...", "sigma=" + scale + " scaled");
	rename("dogTemp1");
	selectImage(baseImageID);
	scale2 = scale * multiplier;
	run("Gaussian Blur...", "sigma=" + scale2 + " scaled");
	rename("dogTemp2");
	imageCalculator("Subtract create 32-bit", "dogTemp1","dogTemp2");
	outputImageID = getImageID;
	rename(outputTitle);
	selectWindow("dogTemp1");
	close();
	selectWindow("dogTemp2");
	close();
	return outputImageID;
}

function addSpotsToOverlay(xc, yc, nSpots, color) {
	// add specified spots to current image's overlay
	run("Select None");
	xc = Array.trim(xc, nSpots);
	yc = Array.trim(yc, nSpots);
	makeSelection("point small " + color + " hybrid", xc, yc);
	Overlay.addSelection();
}

function measureMaskArea(maskID, iROI) {
	// for mask image (0,255) return area in calibrated units
	// use ROI with index 'iROI', or -1 for none
	selectImage(maskID);
	roiManager("select", iROI);
	getRawStatistics(nPixels, mean, min, max, std, hist);
	nForeground = mean * nPixels / 255;
	pixelSize = 1;  // 1 pixel length (to convert to scaled)
	toScaled(pixelSize);  // i.e. pixelSize is now um / pixel
	areaPerPixel = pixelSize * pixelSize;  // 1 pixel area, i.e. um^2
	area_unit2 = nForeground * areaPerPixel;
	run("Select None");
	return area_unit2;
}

function baseName(filename) {
    // return filename string without extension
    if (lastIndexOf(filename, ".") > -1) {
    	return substring(filename, 0, lastIndexOf(filename, "."));
    } else {
    	return filename;
    }
}
