// autoQC_counter2.ijm
// - analysis of outlined cells with 2 markers: classify spots as red or yellow (R+G)
// - 2 separate analyses: point-based (peak detection) and object-based (analyze particles)
//
// Usage:
// - single image and batch modes, selected in dialog
// - both modes expect 2D multi-channel (2+) input image(s)
// - user-specified red-channel (mean +/- N stdDev, for peaks) & ratio thresholds 
// - ROIs for cells *only* should be in ROI manager, or image overlays in batch mode
// - batch mode runs on a folder of .tif images saved from ImageJ with ROI overlays
//
// Output:
// - output images with cell ROIs, multi-point spot ROIs and object ROIs added to overlay
// - summary results: 1 row per cell ROI (results for all images appended to Summary table):
//   - Image name
//   - Cell ROI number 'N', Name and X,Y centroid position
//   - Cell area
//   - number of spots below red threshold intensity, nBlack
//   - number of red spots, nRed
//   - number of red+green spots, nYellow
//   - total areas of: red, yellow, red+yellow
//   - number of red objects, yellow objects
//   - per-cell table with per-particle (object) areas, shape, and 2-channel intensity stats
//   - (parameters: redChannel,greenChannel,rMedianFilter,spotRadius,
//                   minimumObjectArea,peakNstd,rThresh,redMinimumSD,redMinimum)
//
// g.ball@dundee.ac.uk, Dundee Imaging Facility (2023)
// license: MIT License


// Parameters
grnCh = 1;  // channel for grnChondrial marker (green)
redCh = 2;  // channel for mitophagy marker (red)
peakNstd = 2;  // standard deviations above mean for peak detection (red channel)
redMinimumSD = 2.5;  // minimum red intensity for red spots & objects: mean+/-nStdDev
rThresh = 1.0;  // minimum red/green value for red spots
spotRadius = 0.1;  // smallest spot radius, in calibrated units (microns)
minimumObjectArea = 2 * 3.14 * spotRadius * spotRadius;  // minimum area for object (um2)
medianFilterRadius = 1;  // in pixels

// Options
splitTouchingObjects = true;
showRedBlobImage = false;  // option to show red blob image
doBatch = false;  // true = do batch analysis (otherwise, single active image)

// Hard-coded parameters
DoGmultiplier = 5;  // multiplier for difference of gaussian filter
shapeMeasures = "area centroid perimeter fit shape display";
intensMeasures = newArray("Mean", "StdDev", "Min", "RawIntDen", "Median");
amber = "#ffbf00";  // for "yellow" object outlines


// --- start Macro ---

// 1. setup & dialog to check/update channels & options
setOption("BlackBackground", true);  // masks with black background
resetMeasureOptions();
roiManager("Show All with labels");
roiManager("UseNames", "true");  // use ROI name as label
run("Input/Output...", "jpeg=85 gif=-1 file=.csv copy_row copy_column save_column save_row");  // save results as CSV

Dialog.create("autoQC_counter2");
Dialog.addMessage("   Parameters", 16, "#000000");
Dialog.addNumber("Green channel", grnCh);
Dialog.addNumber("Red channel", redCh);
Dialog.addNumber("Median filter radius", medianFilterRadius);
Dialog.addNumber("Spot radius in calibrated units (microns)?", spotRadius);
Dialog.addNumber("Minimum object area (micron^2) ?", minimumObjectArea);
Dialog.addNumber("Number of stdDevs above background for spot detection", peakNstd);
Dialog.addNumber("Minimum red intensity for spots & objects (stdDevs above mean) ", redMinimumSD);
Dialog.addNumber("Minimum red/green intensity ratio for red spots & objects", rThresh);
Dialog.addMessage("   Options", 16, "#000000");
Dialog.addCheckbox("Split touching objects?", splitTouchingObjects);
Dialog.addCheckbox("Show red blob image?", showRedBlobImage);
Dialog.addCheckbox("Batch Mode?", doBatch);
Dialog.show();
grnCh = Dialog.getNumber();
redCh = Dialog.getNumber();
medianFilterRadius = Dialog.getNumber();
spotRadius = Dialog.getNumber();
minimumObjectArea = Dialog.getNumber();
peakNstd = Dialog.getNumber();
redMinimumSD = Dialog.getNumber();
rThresh = Dialog.getNumber();
splitTouchingObjects = Dialog.getCheckbox();
showRedBlobImage = Dialog.getCheckbox();
doBatch = Dialog.getCheckbox();


// 2. for batch mode: get input folder, build list of files, and make output folder
if (doBatch) {
	inputFolder = getDirectory("Choose a folder containing ImageJ .tif images with cell ROI overlay");
	files = getFileList(inputFolder);
	images = selectEnding(files, ".tif");
	numberOfImages = images.length;
	outputFolder = inputFolder + "Results_autoQC_counter2_" + timeStamp();
	print("Saving results in: " + outputFolder);
	File.makeDirectory(outputFolder);
} else {
	numberOfImages = minOf(1, nImages);  // just analyze 1st active image in single-image mode
	if (numberOfImages < 1 || roiManager("count") < 1) {
		exit("Open image & cell ROIs in ROI manager required for single-image mode");
	}
}
// Set up summary results table (all images)
summaryName = "Summary-autoQC2";
Table.create(summaryName);
summaryRow = 0;


// --- 3. Start Iteration over images ---
showProgress(0);
run("Clear Results");
setBatchMode(true);  // hide images during calculations
for (imageNo = 0; imageNo < numberOfImages; imageNo++) {
	
	// 3a. open next image and check
	doAnalyzeThisImage = true;
	if (doBatch) {
		filepath = inputFolder + File.separator + images[imageNo];
		open(filepath);
		imageTitle = getTitle();
		doAnalyzeThisImage = true;
		nCells = Overlay.size;
		if (nCells > 0) {	
			run("To ROI Manager");
		} else {
			print("warning: " + imageTitle + " has no cell ROIs - skipping");
			doAnalyzeThisImage = false;
		}
	} else {
		imageTitle = getTitle();
		nCells = roiManager("count");
		run("Select None");
		run("Duplicate...", "duplicate");  // do not modify original (single image mode)
	}
	Stack.getDimensions(width, height, nc, nz, nt);
	if (nc < 2) {
		if (doBatch) {
			imageTitle = images[imageNo];
		}
		print("warning: " + imageTitle + " has only 1 channel!");
		doAnalyzeThisImage = false;
	}

	// 3b. analyze image
	if (doAnalyzeThisImage) {		
		
		// 3b.I) pre-filter & generate red channel "blob" image
		inpBitDepth = bitDepth();
		Stack.getUnits(xu, yu, zu, tu, vu);
		areaUnits = "" + xu + "2";
		
		// pre-process with median filter
		filteredID = getImageID();
		run("Median...", "radius=1");  // all channels
		filteredTitle = "" + baseName(imageTitle) + "_filtered.tif";
		rename(filteredTitle);
		
		// calculate "redMinimum" threshold red intensity using whole-image mean, SD
		Stack.setChannel(redCh);
		getStatistics(area, redMean, redMin, redMax, redStd);  // whole image
		redMinimum = redMean + (redStd * redMinimumSD);
		
		// generate red blob image (DoG filter) & calculate peakThresh
		redBlobImageID = dogFilter(filteredID, spotRadius, DoGmultiplier);
		rename("redBlobImage");
		selectImage(redBlobImageID);
		run("Select None");
		getStatistics(area, mean, min, max, std);  // whole image
		peakThresh = peakNstd * std;
		
		run("Clear Results");  // default Results table for detailed per-cell per-object results
		
		// 3b.II) iterate over cells - write results row and results image for each
		nRedParticlesPerCell = newArray(nCells);
		nYellowParticlesPerCell = newArray(nCells);
		for (iCell = 0; iCell < nCells; iCell++) {
			cellResultsName = "Results-autoQC2-" + baseName(imageTitle) + "-objects-cell" + iCell;
			selectWindow("redBlobImage");
			Table.set("Image", summaryRow, imageTitle);
			Table.set("RoiN", summaryRow, iCell);
			
			// 3b.II.i) find peaks in DoG-filtered red "spot image"
			roiManager("select", iCell);
			Table.set("RoiName", summaryRow, Roi.getName());
			Table.set("RoiX", summaryRow, getValue("X"));
			Table.set("RoiY", summaryRow, getValue("Y"));
			Table.set("CellArea_" + areaUnits, summaryRow, getValue("Area"));
			run("Duplicate...", "title=redBlobImage-" + iCell);
			run("Find Maxima...", "prominence=" + peakThresh + " output=[Point Selection]");
			getSelectionCoordinates(xc, yc);
			nSpots = xc.length;
			
			// 3b.II.ii) "object" detection & measurement
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
			setThreshold(1, 255);
			areaRedAndYellow = getValue("Area limit");
			resetThreshold();
			
			// calculate mask for above-thresh red spots that are red only (R/G above ratio thresh)
			imageCalculator("AND create", "redAndYellowSpotMask","redOverGreenRatioMask");
			rename("redOnlySpotMask");
			redOnlyMaskID = getImageID();
			
			// threshold & measure red particle area stats
			if (splitTouchingObjects) {
				run("8-bit");  // for watershed
				run("Watershed");
			}
			run("" + inpBitDepth + "-bit");  // for merge
			setThreshold(1, 255);
			areaRedTotal = getValue("Area limit");
			nRoisBefore = roiManager("count");
			run("Set Measurements...", shapeMeasures + " redirect=None decimal=3");
			run("Analyze Particles...", "size=" + minimumObjectArea + "-Infinity display add");
			nRoisAfter = roiManager("count");
			nRedParticlesPerCell[iCell] = nRoisAfter - nRoisBefore;
			redRoi0 = nRoisBefore;
			redRoiN = nRoisAfter-1;
			hadRedObjects = false;
			if (nRedParticlesPerCell[iCell] > 0) {
				hadRedObjects = true;
				rOffset = redRoi0;  // offset between ROI manager indices and Results table rows
				renameRois(redRoi0, redRoiN, "cell" + iCell + "_redObj");
				updateLabelValuesWithROInames(redRoi0, redRoiN, rOffset);
				// measure red object intensity stats
				measureRois("red", redRoi0, redRoiN, greenImageID, intensMeasures, rOffset, "green");
				measureRois("red", redRoi0, redRoiN, redImageID, intensMeasures, rOffset, "red");
			} else {
				rOffset = -1;  // try to update later with first yellow object
			}
			
			// calculate mask image for above-thresh red spots that are "yellow" (R/G below ratio thresh)
			imageCalculator("Subtract create", "redAndYellowSpotMask","redOnlySpotMask");
			rename("yellowSpotMask");
			yellowMaskID = getImageID();
			if (splitTouchingObjects) {
				run("8-bit");  // for watershed
				run("Watershed");
			}
			run("" + inpBitDepth + "-bit");  // for merge
						
			// threshold & measure yellow particle area stats
			setThreshold(1, 255);
			areaYellowTotal = getValue("Area limit");
			nRoisBefore = roiManager("count");
			run("Analyze Particles...", "size=" + minimumObjectArea + "-Infinity display add");
			nRoisAfter = roiManager("count");
			nYellowParticlesPerCell[iCell] = nRoisAfter - nRoisBefore;
			yellowRoi0 = nRoisBefore;
			yellowRoiN = nRoisAfter-1;
			if (nYellowParticlesPerCell[iCell] > 0) {
				renameRois(yellowRoi0, yellowRoiN, "cell" + iCell + "_yellowObj");
				if (rOffset == -1) {
					rOffset = yellowRoi0;  // offset between ROI manager indices and Results table rows (where no red objects)
				}
				updateLabelValuesWithROInames(yellowRoi0, yellowRoiN, rOffset);
				// measure yellow object intensity stats
				measureRois("yellow", yellowRoi0, yellowRoiN, greenImageID, intensMeasures, rOffset, "green");
				measureRois("yellow", yellowRoi0, yellowRoiN, redImageID, intensMeasures, rOffset, "red");
			}
			close("redAndYellowSpotMask");
			
			// rename per-object Results for cell (save & close in batch mode)
			IJ.renameResults(cellResultsName);  // red, yellow object stats for cell
			if (doBatch) {				
				saveAs("Results", outputFolder + File.separator + cellResultsName + ".csv");
				close(cellResultsName + ".csv");
			}
			
			// 3b.II.iii) generate merge image with ROI overlays
			closeImage(ratioMaskID);
			C1name = "c1=redChannelImage";
			C2name = "c2=greenChannelImage";
			C3name = "c3=yellowSpotMask";
			C4name = "c4=redOnlySpotMask";
			channelsString = "" + C1name + " " + C2name + " " +  C3name + " " + C4name + " create";
			run("Merge Channels...", channelsString);
			mergeID = getImageID();
			cellResultsSnapshotTitle = "" + baseName(imageTitle) + "_ROI" + iCell;
			rename(cellResultsSnapshotTitle);
			
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
			addRoiToOverlay(cellRoiIndex, "white", "cell" + iCell);  // cell outline in crop
			addRoisToOverlay(redRoi0, redRoiN, "red");  // red objects
			addRoisToOverlay(yellowRoi0, yellowRoiN, amber);  // yellow objects
			
			// 3b.II.iv) calculate spot stats for this cell ROI
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
			
						
			// 3b.II.v) Add coloured classified spots to overlay
			//          & save snapshot of cell with objects & spot detection overlay
			addSpotsToOverlay(xBlack, yBlack, nBlack, "black", "cell" + iCell + "_blackSpots");
			addSpotsToOverlay(xRed, yRed, nRed, "red", "cell" + iCell + "_redSpots");
			addSpotsToOverlay(xYellow, yYellow, nYellow, "yellow", "cell" + iCell + "_yellowSpots");
			run("Select None");
			autoContrastAllChannels(mergeID);
			if (doBatch) {			
				saveAs("Tiff", outputFolder + File.separator + cellResultsSnapshotTitle + "_objects_and_spots.tif");
				close();  // RGB
			}
			
			
			// 3b.II.vi) write spot stats & parameters to results table row
			Table.set("nBlackSpots", summaryRow, nBlack);
			Table.set("nRedSpots", summaryRow, nRed);
			Table.set("nYellowSpots", summaryRow, nYellow);
			Table.set("areaRed_" + areaUnits, summaryRow, areaRedTotal);
			Table.set("areaYellow_" + areaUnits, summaryRow, areaYellowTotal);
			Table.set("areaRedAndYellow_" + areaUnits, summaryRow, areaRedAndYellow);
			Table.set("nRedObjects", summaryRow, nRedParticlesPerCell[iCell]);
			Table.set("nYellowObjects", summaryRow, nYellowParticlesPerCell[iCell]);
			Table.set("redChannel", summaryRow, redCh);
			Table.set("greenChannel", summaryRow, grnCh);
			Table.set("rMedianFilter", summaryRow, medianFilterRadius);
			Table.set("spotRadius_" + xu, summaryRow, spotRadius);
			Table.set("minimumObjectArea_" + areaUnits, summaryRow, minimumObjectArea);
			Table.set("peakNstd", summaryRow, peakNstd);
			Table.set("peakThresh", summaryRow, peakThresh);
			Table.set("rThresh", summaryRow, rThresh);
			Table.set("redMinimumSD", summaryRow, redMinimumSD);
			Table.set("redMinimum", summaryRow, redMinimum);			
			Table.update();
			summaryRow++;
		}
		
		// 3b.III. clean up before next image / finish
		
		selectImage(redBlobImageID);
		if (showRedBlobImage) {
			run("Select None");
			run("Enhance Contrast", "saturated=0.35");
			if (doBatch) {
				saveAs("Tiff", outputFolder + File.separator + baseName(filteredTitle) + "_red_blob_image.tif");
				close();
			}
		} else {
			close();
		}
	
		if (doBatch) {
			close(filteredTitle);
		} else {
			selectWindow(imageTitle);
			roiManager("Show None");
		}
		
	}
	
	showProgress((imageNo + 1) / numberOfImages);
}

// 4. Save Summary table & exit
if (doBatch) {
	Table.save(outputFolder + File.separator + summaryName + ".csv");
	resetMeasureOptions();
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

function updateLabelValuesWithROInames(iFirst, iLast, offset) {
	// update 'Label' column in Results with ROI labels from roiManager (iFirst..iLast in roiManager)
	// offset is subtracted from ROI index to generate index into results table
	roiManager("deselect");
	for (i = iFirst; i <= iLast; i++) {
		roiManager("select", i);
		setResult("Label", i-offset, Roi.getName());
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

function autoContrastAllChannels(imageID) {
	// auto-contrast all channels for specified image
	selectImage(imageID);
	getDimensions(w, h, nc, nz, nt);
	for (c = 1; c <= nc; c++) {
		Stack.setChannel(c);
	    run("Enhance Contrast", "saturated=0.35");
	}
}

function addSpotsToOverlay(xc, yc, nSpots, color, name) {
	// add specified spots to current image's overlay
	run("Select None");
	xc = Array.trim(xc, nSpots);
	yc = Array.trim(yc, nSpots);
	makeSelection("point small " + color + " hybrid", xc, yc);
	Roi.setName(name);
	Overlay.addSelection();
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

function addRoiToOverlay(roiIndex, color, name) {
	// add a single ROI to active image overlay with specified name and color
	roiManager("deselect");
	roiManager("select", roiIndex);
	Roi.setStrokeColor(color);
	Roi.setName(name);
	run("Add Selection...");
	roiManager("deselect");
}
