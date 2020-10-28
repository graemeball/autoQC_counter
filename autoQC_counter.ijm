// autoQC_counter.ijm
// - analysis of outlined cells with 2 markers: classify spots as red or yellow (R+G)
//
// Usage:
// - single image and batch modes, selected in dialog
// - user-specified red-channel (mean +/- N stdDev, for peaks) & ratio thresholds 
// - both modes expect 2D multi-channel (2+) input image(s)
// - ROIs for cells *only* should be in ROI manager, or image overlays in batch mode
// - batch mode runs on a folder of .tif images saved from ImageJ with ROI overlays
//
// Output:
// - output images with cell ROIs and multi-point spot ROIs added to overlay
// - summary results: 1 row per cell ROI (results for all images appended to table):
//   - Image name
//   - Cell ROI name
//   - ROI XY top-left coords, X_tl & Y_tl
//   - cell area
//   - number of spots below red threshold intensity, nBlack
//   - number of red spots, nRed
//   - number of red+green spots, nYellow
//   - (parameters: peakNstd,redMinimumSD,rThresh,spotSize)
//
// g.ball@dundee.ac.uk, Dundee Imaging Facility (2020)
// license: MIT License


// parameters
doBatch = false;  // true = do batch analysis (otherwise, single active image)
grnCh = 1;  // channel for mitochondrial marker (green)
redCh = 2;  // channel for mitophagy marker (red)
peakNstd = 2;  // standard deviations above mean for peaks (red channel)
redMinimumSD = 1;  // minimum red intensity for spots: mean+/-nStdDev
rThresh = 2;  // minimum red/green value for red spots
spotSize = 0.1;  // calibrated units (microns)
DoGmultiplier = 5;  // multiplier for DoG filter


// --- start Macro ---

// 1. dialog to check/update channels & options
Dialog.create("autoQC_counter");
Dialog.addCheckbox("Batch Mode?", doBatch);
Dialog.addNumber("Green channel", grnCh);
Dialog.addNumber("Red channel", redCh);
Dialog.addNumber("Spot size in calibrated units (microns?)", spotSize);
Dialog.addNumber("Number of stdDevs above background for peaks", peakNstd);
Dialog.addNumber("Minimum red intensity for spots (stdDevs above mean) ", redMinimumSD);
Dialog.addNumber("Minimum red/green intensity ratio for red spots", rThresh);
Dialog.show();
doBatch = Dialog.getCheckbox();
grnCh = Dialog.getNumber();
redCh = Dialog.getNumber();
spotSize = Dialog.getNumber();
peakNstd = Dialog.getNumber();
redMinimumSD = Dialog.getNumber();
rThresh = Dialog.getNumber();


// 2. for batch mode: get input folder, build list of files, and make output folder
if (doBatch) {
	inputFolder = getDirectory("Choose a folder containing ImageJ .tif images with cell ROI overlay");
	files = getFileList(inputFolder);
	images = selectEnding(files, ".tif");
	numberOfImages = images.length;
	outputFolder = inputFolder + "Results_autoQC_counter" + timeStamp();
	print("Saving results in: " + outputFolder);
	File.makeDirectory(outputFolder);
} else {
	numberOfImages = minOf(1, nImages);  // just analyze 1st active image in single-image mode
	if (numberOfImages < 1 || roiManager("count") < 1) {
		exit("Open image & cell ROIs in ROI manager required for single-image mode");
	}
}


// --- 3. Start Iteration over images ---
showProgress(0);
run("Clear Results");
//run("Close All");
resultsRow = 0;  // row in Results table
for (imageNo = 0; imageNo < numberOfImages; imageNo++) {
	setBatchMode(true);

	// 3a. open next image and check
	doAnalyzeThisImage = true;
	if (doBatch) {
		filepath = inputFolder + File.separator + images[imageNo];
		open(filepath);
		inpID = getImageID();
		imageTitle = getTitle();
		doAnalyzeThisImage = true;
		nRois = Overlay.size;
		if (nRois > 0) {	
			run("To ROI Manager");
		} else {
			print("warning: " + imageTitle + " has no cell ROIs - skipping");
			doAnalyzeThisImage = false;
		}
	} else {
		nRois = roiManager("count");
		run("Select None");
		imageTitle = getTitle();
		run("Duplicate...", "duplicate");  // do not modify original stack
		inpID = getImageID();
	}
	Stack.getDimensions(width, height, nc, nz, nt);
	Stack.getUnits(xu, yu, zu, tu, vu);
	if (nc < 2) {
		if (doBatch) {
			imageTitle = images[imageNo];
		}
		print("warning: " + imageTitle + " has only 1 channel!");
		doAnalyzeThisImage = false;
	}

	if (doAnalyzeThisImage) {
		// 3b. pre-process with median filter
		run("Select None");
		run("Median...", "radius=1");  // all channels
	
		// 3c. generate red spot image (DoG filter)
		Stack.setChannel(redCh);
		redSpotImageID = dogFilter(inpID, spotSize, DoGmultiplier);
		rename("redSpotImage");
	
		// 3d. iterate over ROIs (cells) and analyze
		for (r = 0; r < nRois; r++) {
			selectImage(inpID);
			roiManager("select", r);
			Stack.setChannel(redCh);
			// calculate redMinimum using whole-image mean & SD
			getStatistics(area, mean, min, max, std);  // whole image
			redMinimum = mean + (std * redMinimumSD);
			Overlay.addSelection();
			setResult("Image", resultsRow, imageTitle);
			setResult("CellROI", resultsRow, Roi.getName());
			Roi.getBounds(x, y, w, h);
			setResult("X_tl", resultsRow, "" + x);
			setResult("Y_tl", resultsRow, "" + y);
			getStatistics(area, mean, min, max, std);  // cell ROI
			setResult("CellArea_" + xu + "2", resultsRow, area);
			selectImage(redSpotImageID);
			run("Select None");
			getStatistics(area, mean, min, max, std);  // whole image
			peakThresh = peakNstd * std;  
			roiManager("select", r);
			run("Find Maxima...", "prominence=" + peakThresh + " output=[Point Selection]");
			if (selectionType == 10) {
				getSelectionCoordinates(xc, yc);
				nSpots = xc.length;
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
				selectImage(inpID);
				run("Select None");
				iGrn = newArray(nSpots);
				Stack.setChannel(grnCh);
				for (i = 0; i < xc.length; i++) {
					iGrn[i] = getPixel(xc[i], yc[i]);
				}
				iRed = newArray(nSpots);
				Stack.setChannel(redCh);
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
			} else {
				// no spots found
				nBlack = 0;
				nRed = 0;
				nYellow = 0;
			}
			
			// 3e. write stats row to results table
			setResult("nBlack", resultsRow, nBlack);
			setResult("nRed", resultsRow, nRed);
			setResult("nYellow", resultsRow, nYellow);
			setResult("spotSize_" + xu + "2", resultsRow, spotSize);
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

		// 3f. save snapshot of cells with spot detection overlay & clean up
		selectImage(inpID);
		autoContrastAllChannels(inpID);
		if (doBatch) {			
			run("Stack to RGB");
			saveAs("Tiff", outputFolder + File.separator + baseName(images[imageNo]) + "_spots.tif");
			close();  // RGB
		}
		selectImage(redSpotImageID);
		close();
	}
	if (doBatch) {
		selectImage(inpID);
		close();
	}
	showProgress((imageNo + 1) / numberOfImages);
}

// 4. Save Results table & exit
resultsName = "Results-autoQC-redCh" + redCh + "_grnCh" + grnCh + ".csv";
IJ.renameResults(resultsName);
if (doBatch) {
	// configure Results saving: choose .csv
	run("Input/Output...", "jpeg=85 gif=-1 file=.csv copy_row save_column save_row");
	saveAs("Results", outputFolder + File.separator + resultsName);
}

setBatchMode("exit and display");


// --- function definitions ---

function selectEnding(files, ending) {
	// for an array of file names, return a new array of those with specified ending
	withEnding = newArray(0);
	for (i=0; i < files.length; i++) {
		file = files[i];
		if (endsWith(file, ending)) {
			f = newArray(1);
			f[0] = file;
			withEnding = Array.concat(withEnding, f);
		}
	}
	return withEnding;
}

function timeStamp() {
	// generate a time stamp string
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	timeString = toString(year) + "-" + twoDigit(month) + "-" + twoDigit(dayOfMonth);
	DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
	timeString = timeString + "_" + DayNames[dayOfWeek];
	timeString = timeString + twoDigit(hour) + "-" + twoDigit(minute) + "-" + twoDigit(second);
	return timeString;
}

function twoDigit(n) {
	return IJ.pad(n, 2);
}

function baseName(filename) {
	// return filename string without extension
	return substring(filename, 0, lastIndexOf(filename, "."));
}

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

function autoContrastAllChannels(imageID) {
	// auto-contrast all channels for specified image
	selectImage(imageID);
	getDimensions(w, h, nc, nz, nt);
	for (c = 1; c <= nc; c++) {
		Stack.setChannel(c);
	    run("Enhance Contrast", "saturated=0.35");
	}
}
