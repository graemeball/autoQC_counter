# autoQC_counter.ijm

The autoQC counter Fiji/ImageJ macro allows semi-automated batch analysis of autophagy in a set of images. Autophagy is identified using differences in fluorescence intensity between the red mCherry and green GFP proteins of an auto-QC reporter delivered to autophagosomes. 

In batch mode the macro takes as input a folder of images in which an outline for each cell of interest has been created and added to the image overlay in Fiji/ImageJ before saving as tiff. Upon running the macro the user is first asked to select a folder containing images to analyse. A Dialog then prompts the user to check/adjust the following parameters: 1) and 2) green and red channel index, 3) spot size (in microns), 4) a peak finding threshold to identify punctae, 5) a minimum red intensity for punctae of interest, and 6) a red/green intensity ratio threshold to classify punctae as "red" rather than "yellow". Batch analysis of the images in the input folder then begins, and a new output folder is created with a date and timestamp in the name. A new Results table is created for each run, with one row of results per cell region-of-interest, and this is saved to the output folder with a "spots" image corresponding to each input image that shows the location and classification of the punctae identified. Alternatively the macro can be run on the active image only where ROIs for the cells have been added to the ROI manager.

# autoQC_area.ijm
Enhancement to autoQC_counter.ijm for single images only, displaying mask images / areas of "red" vs. "yellow" classification.

