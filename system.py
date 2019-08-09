# export PYTHONPATH=/global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/code/legacypipe/py:$PYTHONPATH
# module load ds9/8.0.1

from legacypipe import queue_calibs
from astropy.io import fits
from contextlib import redirect_stdout
import numpy as np
import os
import sys
import io
        
hdu = fits.open("OutliersList000t035.fits")
#hdu = fits.open("OutliersList042.fits")
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
	ccd = int(input("CCD #: "))
	print("Looking at: ", camera, exposure, ccd)
else:
	index = topTenIndexes[name]
	information = data["Extension name"][index] # the big pixel row with the highest fraction of outliers
	brickname = data["Brick"][index]
	print("Looking at: ", data[index])
	parts = information[0].split('-')
	exposure = int(parts[1])
	ccd = int(parts[2][-1])
	camera = str(parts[0])	
	
os.chdir("/global/project/projectdirs/cosmo/work/legacysurvey/dr8")
hdu = fits.open("survey-ccds-%s-dr8.fits.gz" % (camera))
data = hdu[1].data
ccdS = "CCD" + str(ccd)
counter = 0
surveyLength = hdu[1].header["NAXIS2"]
for elem in range(0, surveyLength):
	try:
		exposureNumber = data[elem]["expnum"]
		ccdNumber = data[elem]["ccdname"]
		if (exposureNumber == exposure) and (ccdNumber == ccdS):
			info = data[elem]
			ccdBRx = data[elem]["crval1"]
			ccdBRy = data[elem]["crval2"]
			ccdCx = data[elem]["ra"]
			ccdCy = data[elem]["dec"]
			ccdX = 2*(ccdCx-ccdBRx)
			ccdY = 2*(ccdCy-ccdBRy)
			#print("CCD Size: ", ccdX, ", ", ccdY) not right?
			filePD = data[elem]["image_filename"]
			info = filePD.split("/")
			date = info[2]
			filename = info[3]
			filenameCut = filename[:-8]
			break
	except IndexError:
		print("Not found!")
		sys.exit()

#print("\npython /global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/code/legacypipe/py/legacypipe/queue-calibs.py --ccds /global/project/projectdirs/cosmo/work/legacysurvey/dr8-garage/zpts/%s/%s/%s-survey.fits --touching" % (camera, date, filenameCut))
#fn = "mosaic/CP20170210/k4m_170211_030552_ooi_zd_v1-survey.fits" # Generate from file

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
for x in range(0, 5):
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
	print(xtl,ytl,xbr,ybr)

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
			raRange = [raL, raU]
			decL = data['DEC1']
			decU = data['DEC2']
			decRange = [decL, decU]
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
	#print(tL[0][4:8])
	if tL[0][4:8] == "p002":
		vertical = 5

tR = str(bR[0][0:4])+str(tL[0][4:8])
bL = str(tL[0][0:4])+str(bR[0][4:8])

print("Top right corner: ", tR)
print("Bottom left corner: ", bL)
print("!!! ", exposure, ccd)

"""
if (vertical == 5):
	if (horizontal == 2) or (horizontal == 3):
		mL2 = str(tL[0][0:5])+"0"+str(np.abs(decTL-2))
		mL3 = str(tL[0][0:5])+"0"+str(np.abs(decTL-3))
		mR2 = str(bR[0][0:5])+"0"+str(np.abs(decTL-2))
		mR3 = str(bR[0][0:5])+"0"+str(np.abs(decTL-3))
		if tL[0][4:8] == "p002":
			mL2 = str(tL[0][0:4])+"p000"
			mR2 = str(bR[0][0:4])+"p000"
		if tL[0][4:8] == "p000":
			mL2 = str(tL[0][0:4])+"m002"
			mR2 = str(bR[0][0:4])+"m002"
		if mL2 in bricks:
			mL = mL2
		else:
			mL = mL3
		if mR2 in bricks:
			mR = mR2
		else:
			mR = mR3
		touchingBricks = [str(bR[0]), mR, tR, bL, mL, str(tL[0])]
		for x in range(0, len(touchingBricks)):
			if len(str(touchingBricks[x])) != 8:
				if len(touchingBricks[x][5:8]) != 3:
					zeroDECa = touchingBricks[x]+"0"
					zeroDECb = touchingBricks[x][0:5]+"0"+touchingBricks[x][-2:]
					#print(zeroDECb)
					if zeroDECb[-4:] == "m003":
						zeroDECb2 = touchingBricks[x][0:4]+"m002"
					if zeroDECb[-4:] == "p003":
						zeroDECb2 = touchingBricks[x][0:4]+"p002"
					if zeroDECa in bricks:
						touchingBricks[x] = zeroDECa
					elif zeroDECb in bricks:
						touchingBricks[x] = zeroDECb
					elif zeroDECb2 in bricks:
						touchingBricks[x] = zeroDECb2
				if len(touchingBricks[x][0:4]) != 4:
					zeroRAb = "0"+touchingBricks[x] 
					zeroRAa = touchingBricks[0:3]+"0"
		print(touchingBricks)
		#print("CCD Covers 6 bricks!")
	if (horizontal == 5):
		tM2 = "0"+str(raTL-2)+str(tL[0][4:8])
		tM3 = "0"+str(raTL-3)+str(tL[0][4:8])
		mL2 = str(tL[0][0:5])+"0"+str(np.abs(decTL-2))
		mL3 = str(tL[0][0:5])+"0"+str(np.abs(decTL-3))
		mM22 = "0"+str(raTL-2)+str(tL[0][4])+"0"+str(np.abs(decTL-2))
		mM23 = "0"+str(raTL-2)+str(tL[0][4])+"0"+str(np.abs(decTL-3))
		mM32 = "0"+str(raTL-3)+str(tL[0][4])+"0"+str(np.abs(decTL-2))
		mM33 = "0"+str(raTL-3)+str(tL[0][4])+"0"+str(np.abs(decTL-3))
		mR2 = str(bR[0][0:5])+"0"+str(np.abs(decTL-2))
		mR3 = str(bR[0][0:5])+"0"+str(np.abs(decTL-3))
		bM2 = "0"+str(raTL-2)+str(bR[0][4:8])
		bM3 = "0"+str(raTL-3)+str(bR[0][4:8])
		if tL[0][4:8] == "p002":
			mM22 = "0"+str(raTL-2)+"p000"
			mM32 = "0"+str(raTL-3)+"p000"
			mR2 = str(bR[0][0:4])+"p000"
			mL2 = str(tL[0][0:4])+"p000"
		if decTL != 3:
			mM22 = mM22[0:5]+"0"+mM22[5:]
			mM23 = mM23[0:5]+"0"+mM23[5:]
			mM32 = mM32[0:5]+"0"+mM32[5:]
			mM33 = mM33[0:5]+"0"+mM33[5:]
		if tM2 in bricks:
			tM = tM2
		else:
			tM = tM3
		if mL2 in bricks:
			mL = mL2
		else:
			mL = mL3
		if mR2 in bricks:
			mR = mR2
		else:
			mR = mR3
		if mM22 in bricks:
			mM = mM22
		elif mM23 in bricks:
			mM = mM23
		elif mM32 in bricks:
			mM = mM32
		elif mM33 in bricks:
			mM = mM33
		if bM2 in bricks:
			bM = bM2
		elif bM3 in bricks:
			bM = bM3
		touchingBricks = [str(bR[0]), mR, tR, bM, mM, tM, bL, mL, str(tL[0])]
		for x in range(0, len(touchingBricks)):
			if len(str(touchingBricks[x])) != 8:
				if len(touchingBricks[x][5:8]) != 3:
					zeroDECa = touchingBricks[x]+"0"
					zeroDECb = touchingBricks[x][0:5]+"0"+touchingBricks[x][-2:]
					if zeroDECa in bricks:
						touchingBricks[x] = zeroDECa
					elif zeroDECb in bricks:
						touchingBricks[x] = zeroDECb
				if len(touchingBricks[x][0:4]) != 4:
					zeroRAb = "0"+touchingBricks[x] 
					zeroRAa = touchingBricks[0:3]+"0"
		#print(touchingBricks)
		#print("CCD Covers 9 Bricks!")
	
elif (vertical == 2) or (vertical == 3):
	if (horizontal == 2) or (horizontal == 3): ## GOOD
		touchingBricks = [str(bR[0]), tR, bL, str(tL[0])]
		#print(touchingBricks)
		#print("CCD Covers 4 bricks!")
	elif horizontal == 5:
		tM2 = "0"+str(raTL-2)+tL[0][4:]
		tM3 = "0"+str(raTL-3)+tL[0][4:]
		bM2 = "0"+str(raTL-2)+bR[0][4:]
		bM3 = "0"+str(raTL-3)+bR[0][4:]
		if raTL != 3:
			tM2 = "0"+tM2
			tM3
		print(bM2, bM3, tM2,tM3)
		if tM2 in bricks:
			tM = tM2
		elif tM3 in bricks:
			tM = tM3
		if bM2 in bricks:
			bM = bM2
		elif bM3 in bricks:
			bM = bM3
		touchingBricks = [str(bR[0]), tR, bM, tM, bL, str(tL[0])]
		print(touchingBricks)
"""
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
#from operator import itemgetter
tmpBricks.sort()

if len(okBricks) == 9:
	touchingBricks = okBricks
	touchingBricks.sort(reverse=True)
	#print("!!!", touchingBricks)
else:
	for b in tmpBricks:
		touchingBricks.append(b[2])

for x in range(0, len(touchingBricks)):
	b3 = touchingBricks[x][0:3]
	dir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr8/north/metrics/%s/" % b3
	os.chdir(dir)
	hdu = fits.open("outlier-mask-%s.fits.fz" % touchingBricks[x])
	#print(dir+"outlier-mask-%s.fits.fz" % touchingBricks[x])
	for elem in range(1, 64):
		try:
			list = hdu[elem].data
			extName = hdu[elem].header['EXTNAME']
			parts = extName.split('-')
			camera = parts[0]
			expnum = int(parts[1])
			ccdNumb = int(parts[2][-1])
			if (expnum == exposure) and (ccdNumb == ccd):	# confirms this brick has same ccd exposure as parsed
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
		os.chdir("/global/project/projectdirs/cosmo/work/legacysurvey/dr8/north/metrics/%s" % b3)
		hdu = fits.open("outlier-mask-%s.fits.fz" % touchingBricks[x])
		#print(touchingBricks[x])
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
	print(exposure, ccd)
	
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
		"""
		print(touchingBricks[5], touchingBricks[3], touchingBricks[1])
		print(touchingBricks[4], touchingBricks[2], touchingBricks[0])
		print(exposure, ccd)
		"""
	
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
	
	"""
	print(touchingBricks[8], touchingBricks[5], touchingBricks[2])
	print(touchingBricks[7], touchingBricks[4], touchingBricks[1])
	print(touchingBricks[6], touchingBricks[3], touchingBricks[0])
	print(exposure, ccd)
	"""
header = fits.Header()
outlierDir = "/global/homes/k/kdevries/outliers"
os.chdir(outlierDir)
fits.writeto('testerYAY.fits', A, header, overwrite=True)
outlierDir = outlierDir+"/testerYAY.fits"
os.system("ds9 "+imdir+" "+outlierDir)