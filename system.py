# Have to run the two lines below in the command line before running system.py
# export PYTHONPATH=/global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/code/legacypipe/py:$PYTHONPATH
# module load ds9/8.0.1

from legacypipe import queue_calibs
from astropy.io import fits
from contextlib import redirect_stdout
import numpy as np
import os
import sys
import io

hdu = fits.open("OutliersList000t035.fits") # OutlierAccess.py returns fits file for certain collection of bricks
data = hdu[1].data
pixels = data['Total pixels']
bigPixels = []
counter = 0
index = []
for x in pixels:
	if x > 100000:
		bigPixels.append(x)
		index.append(counter)
	counter = counter + 1
fractions = data['Fraction of outliers']
fractionBigPixel = []
for i in range(len(index)):
	newF = fractions[index[i]]
	fractionBigPixel.append(newF)

max = np.max(fractionBigPixel)
index = np.where(fractions == max)
max2 = np.max(fractionBigPixel)
maxFractions = fractionBigPixel.copy()

topTenFracs = []
topTenIndexes = []

for frac in range(0, 15):
	max2 = np.max(maxFractions)
	topTenFracs.append(max2)
	ind = np.where(fractions == max2)
	topTenIndexes.append(ind)
	print(str(frac) + ": " + str(data[ind]))
	maxFractions.remove(max2)

print("15: Manual Entry")
name = int(input("Which would you like to choose: "))
if name == 15:
	camera = str(input("Camera: "))
	exposure = int(input("Exposure #: "))
	ccd = input("CCD # (If decam, put in letter and number): ")
	print("Looking at: ", camera, exposure, ccd)
else:
	index = topTenIndexes[name]
	information = data["Extension name"][index] # the big pixel row with the highest fraction of outliers
	brickname = data["Brick"][index]
	print("Looking at: ", data[index])
	parts = information[0].split('-')
	camera = str(parts[0])	
	exposure = int(parts[1])
	if camera == "decam":
		ccd = int(parts[2][1:])
	else:
		ccd = int(parts[2][-1])
	
os.chdir("/global/project/projectdirs/cosmo/work/legacysurvey/dr8")
hdu = fits.open("survey-ccds-%s-dr8.fits.gz" % (camera))
data = hdu[1].data
if camera == "decam":
	ccdS = str(ccd)
	ccdN = ccd[1:]
else:
	ccdS = "CCD" + str(ccd)
	ccdN = ccd
counter = 0
surveyLength = hdu[1].header["NAXIS2"]
for elem in range(0, surveyLength):
	try:
		exposureNumber = data[elem]["expnum"]
		ccdNumber = data[elem]["ccdname"]
		if (exposureNumber == exposure) and (ccdNumber == ccdS):
			info = data[elem]
			filePD = data[elem]["image_filename"]
			info = filePD.split("/")
			date = info[2]
			filename = info[3]
			filenameCut = filename[:-8]
			break
	except IndexError:
		print("Not found!")
		sys.exit()
		
fn = "%s/%s/%s-survey.fits" % (camera, date, filenameCut)
dir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/zpts/%s" % fn

arg1 = "--ccds"
arg2 = "--touching"

f = io.StringIO()
with redirect_stdout(f):
	status = queue_calibs.main([arg1, dir, arg2])
s = f.getvalue()

bricks = []
line = ""
for element in s:
	if element != '\n':
		line += element
	else:
		if len(line) == 8:
			bricks.append(line)
		line = ""
		
ras = []
decs = []
extensions = []
indexes = []
imdir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/zpts/%s" % (fn)
hdu = fits.open(imdir)
data = hdu[1].data
for x in range(0, 100):
	try:
		ccdHere = data[x]['ccdname']
		if str(ccdHere) == str(ccdS):
			ccdElement = x+1
	except IndexError:
		break
imdir = "/global/project/projectdirs/cosmo/staging/%s[%s]" % (filePD, str(ccdElement))
touchingBricks = []

os.chdir("/global/project/projectdirs/cosmo/work/legacysurvey/dr8")
findingSize = fits.open("ccds-annotated-%s-dr8.fits.gz" % camera)
dat = findingSize[1].data
i = np.where((dat['expnum'] == exposure) & (dat['ccdname'] == ccdS))
if camera == "mosaic":
	xtl = dat[i]['ra2']
	ytl = dat[i]['dec2']
	xbr = dat[i]['ra0']
	ybr = dat[i]['dec0']
	
elif camera == "90prime":
	xtl = dat[i]['ra0']
	ytl = dat[i]['dec0']
	xbr = dat[i]['ra2']
	ybr = dat[i]['dec2']

elif camera == "decam":
	xtl = dat[i]['ra1']
	ytl = dat[i]['dec1']
	xbr = dat[i]['ra3']
	ybr = dat[i]['dec3']

hduBricks = fits.open("survey-bricks.fits.gz")
brickData = hduBricks[1].data
for x in range(0, len(bricks)):
	try:	
		if bricks[x] == '':
			break
		else:
			index = np.where(brickData['BRICKNAME'] == bricks[x])
			data = brickData[index]
			raL = data['RA1']
			raU = data['RA2']
			decL = data['DEC1']
			decU = data['DEC2']
			if (xbr < raU) and (xbr > raL):
				if (ybr < decU) and (ybr > decL):
					bR = data['BRICKNAME']
					print("Bottom right corner: ", bR)
			if (xtl < raU) and (xtl > raL):
				if (ytl < decU) and (ytl > decL):
					tL = data['BRICKNAME']
					print("Top left corner: ", tL)
	except IndexError:
			break
					
raTL = int(tL[0][0:4])
decTL = int(tL[0][5:8])
raBR = int(bR[0][0:4])
decBR = int(bR[0][5:8])
if tL[0][4] == "m":
	decTL = -decTL
	decBR = -decBR
	
horizontal = raTL - raBR
vertical = np.abs(decTL-decBR)
if tL[0][4] != bR[0][4]:
	signTL = tL[0][4]
	signBR = bR[0][4]
	if tL[0][4:8] == "p002":
		vertical = 5

tR = str(bR[0][0:4])+str(tL[0][4:8])
bL = str(tL[0][0:4])+str(bR[0][4:8])

#print("Top right corner: ", tR)
#print("Bottom left corner: ", bL)

def brick2coords(brickname):
    ra = 0.1*float(brickname[0:4])
    dec = 0.1*float(brickname[5:8])
    if brickname[4] == "m":
    	dec *= -1.0
    return ra, dec
    
ra1, dec1 = brick2coords(bR[0])
ra2, dec2 = brick2coords(tR)
ra3, dec3 = brick2coords(bL)
ra4, dec4 = brick2coords(tL[0])
# Need to add check for wraparound at ra = 360deg.
tmpBricks = []
okBricks = touchingBricks.copy()
for b in bricks:
	ra, dec = brick2coords(b)
	if ( (ra1 <= ra) and (ra <= ra3) and (dec1 <= dec) and (dec <= dec2) ):
		tmpBricks.append( (ra, dec, b) )
		okBricks.append(b)
tmpBricks.sort()

if len(okBricks) == 9:
	touchingBricks = okBricks
	touchingBricks.sort(reverse=True)
else:
	for b in tmpBricks:
		touchingBricks.append(b[2])
if camera == "decam":
	location = "south"
else:
	location = "north"
for x in range(0, len(touchingBricks)):
	b3 = touchingBricks[x][0:3]
	dir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr8/%s/metrics/%s/" % (location, b3)
	os.chdir(dir)
	hdu = fits.open("outlier-mask-%s.fits.fz" % touchingBricks[x])
	for elem in range(1, 64):
		try:
			list = hdu[elem].data
			extName = hdu[elem].header['EXTNAME']
			parts = extName.split('-')
			camera = parts[0]
			expnum = int(parts[1])
			if camera == "decam":
				ccdNumb = int(parts[2][1:])
			else:
				ccdNumb = int(parts[2][-1])
			if (int(expnum) == int(exposure)) and (int(ccdNumb) == int(ccdN)):	# confirms this brick has same ccd exposure as parsed
				extensions.append(extName)
				indexes.append(elem)
				break
		except IndexError:
			print("Error! The exposure number could not be found for ", touchingBricks[x])
			sys.exit()

dataCorners = []
xs = []
ys = []

for x in range(0, len(touchingBricks)):
	try:
		b3 = touchingBricks[x][0:3]
		os.chdir("/global/project/projectdirs/cosmo/work/legacysurvey/dr8/%s/metrics/%s" % (location, b3))
		hdu = fits.open("outlier-mask-%s.fits.fz" % touchingBricks[x])
		data = hdu[indexes[x]].data
		dataCorners.append(data)
		x0, y0 = data.shape
		xs.append(x0)
		ys.append(y0)
	except IndexError:
		break
		
if len(touchingBricks) == 4:
	nx1 = xs[0] + xs[2]
	nx2 = xs[1] + xs[3]
	ny1 = ys[0] + ys[1]
	ny2 = ys[2] + ys[3]

	if nx1 > nx2:
		nx = nx1
	else:
		nx = nx2
	if ny1 > ny2:
		ny = ny1
	else:
		ny = ny2

	A = np.zeros([nx, ny])
	A[0:xs[0], 0:ys[0]] = dataCorners[0]
	A[0:xs[1], ys[0]:ny1] = dataCorners[1]
	A[xs[0]:nx1, 0:ys[2]] = dataCorners[2]
	A[xs[1]:nx2, ys[2]:ny2] = dataCorners[3]
	print(touchingBricks[3], touchingBricks[1])
	print(touchingBricks[2], touchingBricks[0])
	
elif len(touchingBricks) == 6:
	if vertical == 5:
		nx1 = xs[0] + xs[3]
		nx2 = xs[1] + xs[4]
		nx3 = xs[2] + xs[5]
		ny12 = ys[0] + ys[1] + ys[2]
		ny11 = ys[0] + ys[1]
		ny22 = ys[3] + ys[4] + ys[5]
		ny21 = ys[3] + ys[4]
	
		if nx1 > nx2:
			if nx1 > nx3:
				nx = nx1
			else:
				nx = nx3
		else:
			if nx2 > nx3:
				nx = nx2
			else:
				nx = nx3
	
		if ny12 > ny22:
			ny = ny12
		else:
			ny = ny22
		A = np.zeros([nx,ny])
	
		A[0:xs[0], 0:ys[0]] = dataCorners[0] # Bottom right
		A[0:xs[1], ys[0]:ny11] = dataCorners[1] # Middle right
		A[0:xs[2], ny11:ny12] = dataCorners[2] # Top right
		A[xs[0]:nx1, 0:ys[3]] = dataCorners[3] # Bottom left
		A[xs[1]:nx2, ys[3]:ny21] = dataCorners[4]  # Middle left
		A[xs[2]:nx3, ny21:ny22] = dataCorners[5] # top left
		print(touchingBricks[5], touchingBricks[2])
		print(touchingBricks[4], touchingBricks[1])
		print(touchingBricks[3], touchingBricks[0])
		
	elif (vertical == 2) or (vertical == 3):
		nx11 = xs[0]+xs[2]
		nx12 = xs[0]+xs[2]+xs[4]
		nx21 = xs[1]+xs[3]
		nx22 = xs[1]+xs[3]+xs[5]
		ny1 = ys[0]+ys[1]
		ny2 = ys[2]+ys[3]
		ny3 = ys[4]+ys[5]
		
		if nx12 > nx22:
			nx = nx12
		else:
			nx = nx22
		if ny1 > ny2:
			if ny1 > ny3:
				ny = ny1
			else:
				ny = ny3
		else:
			if ny2 > ny3:
				ny = ny2
			else:
				ny = ny3
		A = np.zeros([nx, ny])
		A[0:xs[0], 0:ys[0]] = dataCorners[0]
		A[0:xs[1], ys[0]:ny1] = dataCorners[1]
		A[xs[0]:nx11, 0:ys[2]] = dataCorners[2]
		A[xs[1]:nx21, ys[2]:ny2] = dataCorners[3]
		A[nx11:nx12, 0:ys[4]] = dataCorners[4]
		A[nx21:nx22, ys[4]:ny3] = dataCorners[5]

		print(touchingBricks[5], touchingBricks[3], touchingBricks[1])
		print(touchingBricks[4], touchingBricks[2], touchingBricks[0])
	
elif len(touchingBricks) == 9:
	nx11 = xs[0]+xs[3]
	nx21 = xs[1]+xs[4]
	nx31 = xs[2]+xs[5]
	nx12 = xs[0]+xs[3]+xs[6]
	nx22 = xs[1]+xs[4]+xs[7]
	nx32 = xs[2]+xs[5]+xs[8]
	ny11 = ys[0]+ys[1]
	ny21 = ys[3]+ys[4]
	ny31 = ys[6]+ys[7]
	ny12 = ys[0]+ys[1]+ys[2]
	ny22 = ys[3]+ys[4]+ys[5]
	ny32 = ys[6]+ys[7]+ys[8]
	
	if nx12 > nx22:
		if nx12 > nx32:
			nx = nx12
		else:
			nx = nx32
	else:
		if nx22 > nx32:
			nx = nx22
		else:
			nx = nx32
	if ny12 > ny22:
		if ny12 > ny32:
			ny = ny12
		else:
			ny = ny32
	else:
		if ny22 > ny32:
			ny = ny22
		else:
			ny = ny32
			
	A = np.zeros([nx,ny])
	
	A[0:xs[0], 0:ys[0]] = dataCorners[0] # Bottom right
	A[0:xs[1], ys[0]:ny11] = dataCorners[1] # Middle right
	A[0:xs[2], ny11:ny12] = dataCorners[2] # Top right
	A[xs[0]:nx11, 0:ys[3]] = dataCorners[3] # bottom middle
	A[xs[1]:nx21, ys[3]:ny21] = dataCorners[4] # middle middle
	A[xs[2]:nx31, ny21:ny22] = dataCorners[5] # top middle
	A[nx11:nx12, 0:ys[6]] = dataCorners[6] # bottom left
	A[nx21:nx22, ys[6]:ny31] = dataCorners[7] # middle left
	A[nx31:nx32, ny31:ny32] = dataCorners[8] # top left
	
	print(touchingBricks[0], touchingBricks[3], touchingBricks[6])
	print(touchingBricks[1], touchingBricks[4], touchingBricks[7])
	print(touchingBricks[2], touchingBricks[5], touchingBricks[8])

header = fits.Header()
outlierDir = "/global/homes/k/kdevries/outliers"
os.chdir(outlierDir)
fits.writeto('testerYAY.fits', A, header, overwrite=True)
outlierDir = outlierDir+"/testerYAY.fits"
os.system("ds9 "+imdir+" "+outlierDir)