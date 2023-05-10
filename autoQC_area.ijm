// autoQC_area.ijm
// - enhancement to autoQC_counter.ijm for single images only, displaying area masks 
//
// Usage:
// - user-specified red-channel (mean +/- N stdDev, for peaks) & ratio thresholds 
// - both modes expect 2D multi-channel (2+) input image(s)
// - ROIs for cells *only* should be in ROI manager or image overlay
//
// Output:
// - output image with cell ROIs and multi-point spot ROIs added to overlay:
//   - with additional "yellow object" and "red object" output channels
//   - output channel 1=red, 2=green, 3=yellow objects (yellow), 4=red objects (magenta)
// - summary results: 1 row per cell ROI (results for all images appended to table):
//   - Image name
//   - Cell ROI name and X,Y centroid position
//   - Cell area
//   - number of spots below red threshold intensity, nBlack
//   - number of red spots, nRed
//   - number of red+green spots, nYellow
//   - total areas of: red, yellow, red+yellow
//   - per-cell table with per-particle areas, shape, and 2-channel intensity stats
//   - (parameters: peakNstd,redMinimumSD,rThresh,spotSize)
//
// g.ball@dundee.ac.uk, Dundee Imaging Facility (2020-2023)
// license: MIT License


// parameters
grnCh = 1;  // channel for grnChondrial marker (green)
redCh = 2;  // channel for mitophagy marker (red)
peakNstd = 2;  // standard deviations above mean for peaks (red channel)
redMinimumSD = 2.5;  // minimum red intensity for spots: mean+/-nStdDev
rThresh = 1.0;  // minimum red/green value for red spots
spotSize = 0.1;  // smallest spot radius, in calibrated units (microns)
DoGmultiplier = 5;  // multiplier for DoG filter
showRedBlobImage = false;  // option to show red blob image
shapeMeasures = "area centroid perimeter fit shape";
intensMeasures = newArray("Mean", "StdDev", "Min", "RawIntDen", "Median");
amber = "#ffbf00";


// --- start Macro ---
if (Overlay.size > 0) {
	run("To ROI Manager");
}

// check we have at least 1 ROI & composite image before going further
nCells = roiManager("count");
if (nCells < 1) {
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

// image summary results table
paramString = "R" + redCh + "_G" + grnCh;
resultsName = "Results-autoQC-" + paramString;
Table.create(resultsName);

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
resetMeasureOptions();

// pre-process with median filter
run("Select None");
run("Duplicate...", "duplicate");
filteredID = getImageID();
run("Median...", "radius=1");  // all channels
filteredTitle = "" + baseName(inpTitle) + "_filtered.tif";
rename(filteredTitle);

// calculate "redMinimum" threshold red intensity using whole-image mean, SD
Stack.setChannel(redCh);
getStatistics(area, redMean, redMin, redMax, redStd);  // whole image
redMinimum = redMean + (redStd * redMinimumSD);

// generate red blob image (DoG filter) & calculate peakThresh
redBlobImageID = dogFilter(filteredID, spotSize, DoGmultiplier);
rename("redBlobImage");
selectImage(redBlobImageID);
run("Select None");
getStatistics(area, mean, min, max, std);  // whole image
peakThresh = peakNstd * std;
minimumSpotVolume = 2 * 3.14 * spotSize * spotSize;  // >smallest volume (um2)

run("Clear Results");
resultsRow = 0;
// iterate over cells - write results row and results image for each
nRedParticlesPerCell = newArray(nCells);
nYellowParticlesPerCell = newArray(nCells);
for (iCell = 0; iCell < nCells; iCell++) {
	cellResultsName = resultsName + "-cell" + iCell;
	selectWindow("redBlobImage");
	Table.set("Image", resultsRow, inpTitle);
	Table.set("RoiN", resultsRow, iCell);
	// find peaks in DoG-filtered red "spot image"
	roiManager("select", iCell);
	Table.set("RoiName", resultsRow, Roi.getName());
	Table.set("RoiX", resultsRow, getValue("X"));
	Table.set("RoiY", resultsRow, getValue("Y"));
	Table.set("CellArea_" + areaUnits, resultsRow, getValue("Area"));
	run("Duplicate...", "title=redBlobImage-" + iCell);
	run("Find Maxima...", "prominence=" + peakThresh + " output=[Point Selection]");
	getSelectionCoordinates(xc, yc);
	nSpots = xc.length;
	
	// create red spot image mask (before intensity & ratio thresholds)
	roiManager("select", iCell);
	setThreshold(peakThresh, max);
	run("Convert to Mask");
	rename("redSpotMask");
	redSpotMaskID = getImageID();
	
	// create red/green ratio image mask
	selectImage(filteredID);
	roiManager("select", iCell);
	Stack.setChannel(redCh);
	run("Duplicate...", " ");
	roiManager("Add");
	cellRoiIndex = roiManager("count") - 1;
	
	rename("redChannelImage");
	redImageID = getImageID();
	selectImage(filteredID);
	roiManager("select", iCell);
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
	
	// create redMinimumMask
	selectImage(redImageID);
	run("Duplicate...", " ");
	redMax = getValue("Max");
	setThreshold(redMinimum, redMax);
	run("Convert to Mask");
	rename("redAboveMinimumMask");
	redAboveMinimumID = getImageID();
	
	// calculate mask for above-threshold red spots (i.e. red & yellow spots)
	imageCalculator("AND create", "redSpotMask","redAboveMinimumMask");
	rename("redAndYellowSpotMask");
	redAndYellowSpotMaskID = getImageID();
	closeImage(redSpotMaskID);
	closeImage(redAboveMinimumID);
	selectImage(redAndYellowSpotMaskID);
	areaRedAndYellow = getValue("Area");
	
	// calculate mask for above-thresh red spots that are red only (R/G above ratio thresh)
	imageCalculator("AND create", "redAndYellowSpotMask","redOverGreenRatioMask");
	rename("redOnlySpotMask");
	redOnlyMaskID = getImageID();
	
	// threshold & measure red particle area stats
	run("8-bit");  // for watershed
	run("Watershed");
	run("" + inpBitDepth + "-bit");  // for merge
	setThreshold(1, 255);
	areaRedTotal = getValue("Area limit");
	nRoisBefore = roiManager("count");
	run("Set Measurements...", shapeMeasures + " redirect=None decimal=3");
	run("Analyze Particles...", "size=" + minimumSpotVolume + "-Infinity display add");
	nRoisAfter = roiManager("count");
	nRedParticlesPerCell[iCell] = nRoisAfter - nRoisBefore;
	redRoi0 = nRoisBefore;
	redRoiN = nRoisAfter-1;
	hadRedObjects = false;
	if (redRoiN - redRoi0 > 0) {
		hadRedObjects = true;
		rOffset = redRoi0;  // offset between ROI manager indices and results table rows
		renameRois(redRoi0, redRoiN, "cell" + iCell + "_redObj");
		// measure red object intensity stats
		measureRois("red", redRoi0, redRoiN, greenImageID, intensMeasures, rOffset, "green");
		measureRois("red", redRoi0, redRoiN, redImageID, intensMeasures, rOffset, "red");
	}
	
	// calculate mask image for above-thresh red spots that are "yellow" (R/G below ratio thresh)
	imageCalculator("Subtract create", "redAndYellowSpotMask","redOnlySpotMask");
	rename("yellowSpotMask");
	yellowMaskID = getImageID();
	run("8-bit");  // for watershed
	run("Watershed");
	run("" + inpBitDepth + "-bit");  // for merge
	
	// threshold & measure yellow particle area stats
	setThreshold(1, 255);
	areaYellowTotal = getValue("Area limit");
	nRoisBefore = roiManager("count");
	run("Analyze Particles...", "size=" + minimumSpotVolume + "-Infinity display add");
	nRoisAfter = roiManager("count");
	nYellowParticlesPerCell[iCell] = nRoisAfter - nRoisBefore;
	yellowRoi0 = nRoisBefore;
	yellowRoiN = nRoisAfter-1;
	if (yellowRoiN - yellowRoi0 > 0) {
		renameRois(yellowRoi0, yellowRoiN, "cell" + iCell + "_yellowObj");
		// measure yellow object intensity stats
		measureRois("yellow", yellowRoi0, yellowRoiN, greenImageID, intensMeasures, rOffset, "green");
		measureRois("yellow", yellowRoi0, yellowRoiN, redImageID, intensMeasures, rOffset, "red");
	}
	close("redAndYellowSpotMask");
	
	
	IJ.renameResults(cellResultsName);  // red, yellow object stats for cell
	
	// generate merge image with ROI overlays
	run("" + inpBitDepth + "-bit");  // convert to input bit-depth for merge
	closeImage(ratioMaskID);
	//selectImage(greenImageID);
	C1name = "c1=redChannelImage";
	C2name = "c2=greenChannelImage";
	//C3name = "c3=redAndYellowSpotMask";
	C3name = "c3=yellowSpotMask";
	C4name = "c4=redOnlySpotMask";
	channelsString = "" + C1name + " " + C2name + " " +  C3name + " " + C4name + " create";
	run("Merge Channels...", channelsString);
	mergeID = getImageID();
	rename(filteredTitle + "_ROI" + iCell);
	
	Stack.setChannel(1);
	run("Red");
	appendChannelInfoToSliceLabels(mergeID, 1, "red");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(2);
	run("Green");
	appendChannelInfoToSliceLabels(mergeID, 2, "green");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(3);
	run("Yellow");
	appendChannelInfoToSliceLabels(mergeID, 3, "yellowObjects");
	setMinAndMax(0, 255);
	Stack.setChannel(4);
	run("Magenta");
	appendChannelInfoToSliceLabels(mergeID, 4, "redObjects");
	setMinAndMax(0, 255);
	addRoisToOverlay(cellRoiIndex, cellRoiIndex, "white");  // cell outline in crop
	addRoisToOverlay(redRoi0, redRoiN, "red");  // red objects
	addRoisToOverlay(yellowRoi0, yellowRoiN, amber);  // yellow objects
	// calculate spot stats for this cell ROI
	selectImage(mergeID);
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
	// Add coloured classified spots to overlay
	addSpotsToOverlay(xBlack, yBlack, nBlack, "black");
	addSpotsToOverlay(xRed, yRed, nRed, "red");
	addSpotsToOverlay(xYellow, yYellow, nYellow, "yellow");
	run("Select None");
	
	// write stats row to results table
	Table.set("nBlackSpots", resultsRow, nBlack);
	Table.set("nRedSpots", resultsRow, nRed);
	Table.set("nYellowSpots", resultsRow, nYellow);
	Table.set("areaRed_" + areaUnits, resultsRow, areaRedTotal);
	Table.set("areaYellow_" + areaUnits, resultsRow, areaYellowTotal);
	Table.set("areaRedAndYellow_" + areaUnits, resultsRow, areaRedAndYellow);
	Table.set("nRedObjects", resultsRow, nRedParticlesPerCell[iCell]);
	Table.set("nYellowObjects", resultsRow, nYellowParticlesPerCell[iCell]);
	Table.set("spotSize_" + areaUnits, resultsRow, spotSize);
	Table.set("peakNstd", resultsRow, peakNstd);
	Table.set("peakThresh", resultsRow, peakThresh);
	Table.set("rThresh", resultsRow, rThresh);
	Table.set("redMinimumSD", resultsRow, redMinimumSD);
	Table.set("redMinimum", resultsRow, redMinimum);
	updateResults();
	resultsRow++;
}

close(filteredTitle);
selectImage(redBlobImageID);
if (showRedBlobImage) {
	run("Select None");
	run("Enhance Contrast", "saturated=0.35");
} else {
	close();
}

setBatchMode("exit and display");

selectWindow(inpTitle);
roiManager("Show None");


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

function appendChannelInfoToSliceLabels(imageID, channel, channelInfo) {
    // for given hyperstack channel, append info to end of each slice label
    selectImage(imageID);
    getDimensions(w, h, nc, nz, nt);
    //setBatchMode("hide");
    for (z = 1; z <= nz; z++) {
        for (t = 1; t <= nt; t++) {
            Stack.setPosition(channel, z, t); 
            sliceLabel = Property.getSliceLabel();
            Property.setSliceLabel(sliceLabel + channelInfo);
        }   
    }   
    //setBatchMode("exit and display");
}

//function id2ti(imageID) {
//	// return image title for given imageID
//	startingID = getImageID();
//	selectImage(imageID);
//	imageTitle = getTitle();
//	selectImage(startingID);
//	return imageTitle;
//}

function closeImage(imageID) {
	// close image with given imageID
	selectImage(imageID);
	close();
}

function resetMeasureOptions() {
	// reset measurement options to defaults (N.B. this blanks other measurements!)
	defaultOptions = "area mean standard min max centroid integrated display";
	run("Set Measurements...", defaultOptions + " redirect=None decimal=3");	
}

function renameRois(iFirst, iLast, prefix) {
	// rename a list of ROIs from iFirst to iLast inclusive, to: prefix0..prefixN
	roiManager("deselect");
	n = 0;
	for (i = iFirst; i <= iLast; i++) {
		roiManager("select", i);
		roiManager("rename", prefix + n);
		n++;
	}
	roiManager("deselect");
}

function addRoisToOverlay(iFirst, iLast, color) {
	// add ROIs iFirst to iLast inclusive to Overlay with color specified
	roiManager("deselect");
	for (i = iFirst; i <= iLast; i++) {
		roiManager("select", i);
		Roi.setStrokeColor(color);
		run("Add Selection...");
	}
	roiManager("deselect");
}

function measureRois(object, iFirst, iLast, imageID, measurements, offset, label) {
	// measure series of ROIs on current channel of active image, where:
	//  object added as new column to identify what ROIs are
	//  iFirst, iLast are first and last ROI indices of interest in manager
	//  imageID is the ID of the image to measure
	//  measurements is list of measurements for getValue()
	//  offset is subtracted from ROI index to generate index into results table
	//  label is prefixed to measurement column names to indicate channel measured
	//  results are written to default Results table 
	selectImage(imageID);
	for (i = iFirst; i <= iLast; i++) {
		roiManager("select", i);
		setResult("Object", i-offset, object);
		for (m = 0; m < measurements.length; m++) {
			setResult(label+measurements[m], i-offset, getValue(measurements[m]));
		}
	}
}
