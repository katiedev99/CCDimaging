# cd /global/project/projectdirs/cosmo/work/legacysurvey/dr8/north/metrics/016
# scp OutlierAccess.py kdevries@cori.nersc.gov:/global/homes/k/kdevries/outliers

from astropy.io import fits
import numpy as np
import glob

workdir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr8/north/metrics/"
brickName = []
exposureMax = []
extnameArray = []
extNArray = []
totPixelsArr = []
negOutliers = []
posOutliers = []
fractionOut = []
totalExtensions = 0
totalImages = 0
totalFolders = 0

# grab command line argument
import argparse
#import sys
parser = argparse.ArgumentParser()
parser.add_argument("sd", type=str, help="First three digits of subdirectory")
args = parser.parse_args()
#print(args.sd)
#sys.exit()

for name in glob.glob(workdir+args.sd+'/outlier-mask-*.fits.fz'):
#for name in glob.glob(workdir+'/outlier-mask-*.fits.fz'):
	hdu = fits.open(name)
	extensionNumber = len(hdu)
	brick = hdu[0].header['BRICK']
	extensionN = 0
	totalImages = totalImages + 1

	for x in range(1, extensionNumber):
		extensionN = extensionN + 1
		extNArray.append(extensionN)
		extname = hdu[x].header['EXTNAME'] 
		extnameArray.append(extname)
		brickName.append(brick)

		pixel1 = hdu[x].header['NAXIS1']
		pixel2 = hdu[x].header['NAXIS2']
		totalPixels = pixel1*pixel2
		totPixelsArr.append(totalPixels)

		data = hdu[x].data
		np.unique(data, return_counts=True)
		unique, counts = np.unique(data, return_counts = True)
		dict(zip(unique, counts))
		zeros = counts[0]
		ones = 0
		twos = 0
		fracOutliers = 0
		if len(unique) >= 2:
			if unique[1] == 1:
				ones = counts[1]
			else:
				twos = counts[1]
			if len(unique) == 3:
				ones = counts[1]
				twos = counts[2]
		fracOutliers = (ones + twos) / totalPixels
		posOutliers.append(ones)
		negOutliers.append(twos)
		fractionOut.append(fracOutliers)
		totalExtensions = totalExtensions+1
		if (totalExtensions % 10000 == 0):
			print("Extensions parsed: " , totalExtensions)
	maximum = np.amax(fractionOut)
	if maximum == 0.0:
		i = None
	else:
		i = fractionOut.index(maximum) + 1
	if (totalImages % 50 == 0):
		totalFolders = totalFolders + 1
		print("Total folders parsed: ", totalFolders)

c1 = fits.Column(name= "Brick", array=brickName, format="8A")
c2 = fits.Column(name= "Extension Name", array=extnameArray, format="35A")
c3 = fits.Column(name= "Extension number", array=extNArray, format="K")
c4 = fits.Column(name= "Total pixels", array=totPixelsArr, format="K")
c5 = fits.Column(name= "Positive outliers", array=posOutliers, format="K")
c6 = fits.Column(name= "Negative outliers", array=negOutliers, format="K")
c7 = fits.Column(name= "Fraction of outliers", array=fractionOut, format="E")
table = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7])
table.writeto('OutliersList'+args.sd+'.fits', overwrite=True)
print("Total number of extensions collected: ", totalExtensions)
print("Total number of images parsed: ", totalImages)