from itertools import product
import os, re
import json

from fiji.util.gui import GenericDialogPlus as GenericDialog
from ij.io import Opener
import ij
from ij import ImagePlus, ImageStack, IJ, WindowManager
from ij.plugin import HyperStackConverter, ChannelSplitter, ImageCalculator
from ij.plugin.filter import BackgroundSubtracter, Analyzer
from ij.plugin.frame import RoiManager
from ij.process import StackStatistics, AutoThresholder
from ij.measure import ResultsTable


methodsList = ["Default",
"Huang",
"IJ_IsoData",
"Intermodes",
"IsoData" ,
"Li" ,
"MaxEntropy",
"Mean",
"MinError",
"Minimum",
"Moments",
"Otsu",
"Percentile",
"RenyiEntropy",
"Shanbhag",
"Triangle",
"Yen"]

def getCombos(dicts, groups):
	retval = {}
	for d in dicts:
	  idTuple = tuple(d[g] for g in groups)
	  z = retval.get(idTuple)
	  if z is None: retval[idTuple] = []
	  retval[idTuple].append(d)
	return retval


#This function is necessary because Fiji's StackStatistics.median is buggy and always returns 0
def calcMedian(ip):
  hist = StackStatistics(ip).histogram16
  histsum = sum(hist)
  ii = 0
  for jj in range(0, 65536):
    ii += hist[jj]
    if ii > histsum / 2:
      break
  return jj


#This section defines image opening and threshold methods
def openUnstitched(group):
  print(group)
  ipArr = []
  frameNo = 1
  frameNoLink = {}
  nucIpArr = []
  cmIpArr = []
  for g in group:
    checkKeys = {k:g[k] for k in g if k not in ["channel", "fileName"]}
    existsNuc = False
    existsCm = False
    for gg in group:
  	  checkKeys2 = {k:g[k] for k in g if k not in ["channel", "fileName"]}
  	  if checkKeys2:
  	    if gg['channel'] == nucChannel: existsNuc = True
  	    elif gg['channel'] == cmChannel: existsCm = True
    if not (existsCm and existsNuc): continue
    fileName = g['fileName']
    channel = g['channel']
    checkKeys
    if channel == nucChannel: nucIpArr.append(opener.openImage(folderPath + fileName))
    elif channel == cmChannel: cmIpArr.append(opener.openImage(folderPath + fileName))
  nucStack = ImageStack(nucIpArr[0].width, cmIpArr[0].height)
  cmStack = ImageStack(nucIpArr[0].width, cmIpArr[0].height)
  for nucIp, cmIp in nucIpArr, cmIpArr:
    nucStack.addSlice(nucIp.getTitle(), nucIp.getProcessor())
    cmStack.addSlice(cmIp.getTitle(), cmIp.getProcessor())
  nucIp = ImagePlus("images{}".format(fileName), nucIpArr[0].getProcessor())
  cmIp = ImagePlus("images{}".format(fileName), cmIpArr[0].getProcessor())
  nucIp.setStack(nucStack)
  cmIp.setStack(cmStack)
  return nucIp, cmIp

def openStitched(group, rowNo, colNo):
  nucch = ""
  cellch = ""
  for g in group:
    if g['channel'] == nucChannel: nucch = g['fileName']
    elif g['channel'] == cmChannel: cellch = g['fileName']
  if nucch == "" or cellch == "":
    print("Match not Found", group)
    return None, None
  nucIp = opener.openImage(folderPath + nucch)
  cmIp = opener.openImage(folderPath + cellch)
  IJ.run(nucIp, "Montage to Stack...", "columns={} rows={} border=0".format(colNo, rowNo))
  nucIp.close()
  nucIp = WindowManager.getImage("Stack")
  nucIp.setTitle("NuclearStain")
  cmIp = IJ.run(cmIp, "Montage to Stack...", "columns={} rows={} border=0".format(colNo, rowNo))
  cmIp = WindowManager.getImage("Stack")
  cmIp.setTitle("CmStain")
  nucIp.show(), cmIp.show()
  return nucIp, cmIp


def generateNucleusMask(nucIp, method="Otsu", stack=True, minNucleusSize=50):
  nucMaskIp = nucIp.duplicate()
  if not stack:
  	  IJ.setAutoThreshold(nucMaskIp, "{} dark".format(method))
	  IJ.run(nucMaskIp, "Convert to Mask", "method={} background=Dark calculate".format(method))
  else:
  	  IJ.setAutoThreshold(nucMaskIp, "{} dark stack".format(method))
	  IJ.run(nucMaskIp, "Convert to Mask", "method={} background=Dark stack".format(method))
  IJ.run(nucMaskIp, "Open", "stack")
  IJ.run("Set Measurements...", "stack redirect=None decimal=3")
  IJ.run(nucMaskIp, "Analyze Particles...", "size={}-Infinity circularity=0.2-1.00 show=Masks exclude in_situ stack".format(minNucleusSize))
  IJ.run(nucMaskIp,"Fill Holes", "stack")
  nucMaskIp.show()
  return nucMaskIp

def generateCardiomyocyteMask(cmIp, method="Otsu", stack=True, minCmSize=500, brightfield = False):
  cmIp.setTitle("CmStain")
  if brightfield:
    IJ.run(cmIp, "Find Edges", "stack")
  if not stack:
    IJ.setAutoThreshold(cmIp, "{} dark".format(method))
    IJ.run(cmIp, "Convert to Mask", "method={} background=Dark calculate".format(method))
  else:
    IJ.setAutoThreshold(cmIp, "{} dark stack".format(method))
    IJ.run(cmIp, "Convert to Mask", "method={} background=Dark stack".format(method))
  if brightfield:
    IJ.run(cmIp, "Morphological Filters (3D)", "operation=Closing element=Ball x-radius=2 y-radius=2 z-radius=0")
    delIp = cmIp
    delIp.changes = False
    delIp.close()
    cmIp = WindowManager.getImage("CmStain-Closing")
    IJ.run(cmIp,"Fill Holes", "stack")
  cmIp.show()
  IJ.run("Set Measurements...", "stack redirect=None decimal=3")
  IJ.run(cmIp, "Analyze Particles...", "size={}-Infinity show=[Masks] exclude in_situ stack".format(cmMinSize))
  cmMaskIp = cmIp.duplicate()
  IJ.run(cmIp, "Analyze Particles...", "size={}-Infinity show=[Count Masks] exclude in_situ stack".format(cmMinSize))
  cmMaskIp.show(), cmIp.show()
  return cmIp, cmMaskIp

import re
preGd = GenericDialog("Output Folder")
preGd.addDirectoryField("Output folder", "")
preGd.showDialog()
outputFolder = preGd.getNextString()
os.chdir(outputFolder)
folderPath, formatString, groupBy, nucChannel, cmChannel, stitched, analyzeNucStack, brightfield, nucMethod, cmMethod, analyzeCmStack, rowNo, colNo, nucMinSize, cmMinSize =\
"",          "",          "",      "DAPI",      "GFP",    False,    True,            False,      "Otsu",    "Triangle", True,         8,     7,      50,        500
__dict__ = globals()
try:
	with open("savedSettings.json") as f:
	  thing = json.load(f)
	  for k, v in thing.iteritems():
	    __dict__[k] = v
	print(thing)
except Exception as e: print(e)
gd = GenericDialog("Analysis parameters")
gd.addDirectoryField("Image folder location", folderPath )
gd.addStringField("Image filename pattern", formatString, 70)
gd.addStringField("Group By", groupBy, 70)
gd.addStringField("Nuclear Stain Channel Name", nucChannel, 70)
gd.addStringField("Cardiomyocyte Image Channel Name", cmChannel, 70)
gd.addCheckbox("Are images stitched by well", stitched)
gd.addChoice("Nuclear thresholding method", methodsList, nucMethod)
gd.addCheckbox("Threshold nuclei by stack?", analyzeNucStack)
gd.addCheckbox("Are the cardiomyocyte images brightfield?", brightfield)
gd.addChoice("Cardiomyocyte thresholding method", methodsList,  cmMethod)
gd.addCheckbox("Threshold cardiomyocyes by stack?", analyzeCmStack)
gd.addNumericField("Rows per well", rowNo, 0)
gd.addNumericField("Columns per well", colNo, 0)
gd.addNumericField("Nuclear minimum size (pixels)", nucMinSize, 0)
gd.addNumericField("Cardiomyocyte minimum size (pixels)", cmMinSize, 0)
gd.showDialog()
folderPath = gd.getNextString() + "/"
formatString = gd.getNextString()
groupBy = re.findall(r"\w+", gd.getNextString())
nucChannel = gd.getNextString()
cmChannel = gd.getNextString()
stitched = gd.getNextBoolean()
nucMethod = gd.getNextChoice()
analyzeNucStack = gd.getNextBoolean()
brightfield = gd.getNextBoolean()
cmMethod = gd.getNextChoice()
analyzeCmStack = gd.getNextBoolean()
rowNo = int(gd.getNextNumber())
colNo = int(gd.getNextNumber())
nucMinSize = int(gd.getNextNumber())
cmMinSize = int(gd.getNextNumber())
allFileNames = os.listdir(folderPath)
jsonStoreDict = {"folderPath": folderPath,
				"formatString":formatString,
				"groupBy": ','.join(groupBy),
				"nucChannel": nucChannel,
				"cmChannel": cmChannel,
				"stitched": stitched,
				"nucMethod": nucMethod,
				"analyzeNucStack": analyzeNucStack,
				"brightfield": brightfield,
				"cmMethod": cmMethod,
				"analyzeCmStack": analyzeCmStack,
				"rowNo": rowNo,
				"colNo": colNo,
				"nucMinSize": nucMinSize,
				"cmMinSize": cmMinSize,
				"allFileNames": allFileNames}
with open("savedSettings.json", "w+") as f:
  json.dump(jsonStoreDict, f)
reList = []
for fileName in allFileNames:
  match = re.match(formatString, fileName)
  if match is None: continue
  matchDict = match.groupdict()
  matchDict['fileName'] = match.string
  reList.append(matchDict)
opener = Opener()
bs = BackgroundSubtracter()
for outerPairs, group in getCombos(reList, groupBy).iteritems():
  print(group)
  nucIp, cmIp = openStitched(group, rowNo, colNo) if stitched else openUnstitched(group)
  if nucIp is None or cmIp is None: continue
  nucStack, cmStack = nucIp.getStack(), cmIp.getStack()
  nucIpMedian = calcMedian(nucIp)
  IJ.run(nucIp, "Subtract...", "value=" + str(nucIpMedian) + " stack")
  nucIpForMeasure = nucIp.duplicate()
  nucMaskIp = generateNucleusMask(nucIp, nucMethod, analyzeNucStack, nucMinSize)
  cmIp, cmMaskIp = generateCardiomyocyteMask(cmIp, cmMethod, analyzeCmStack, cmMinSize, brightfield)
  nucMaskIp.show()
  nucIpForMeasure.show()
  rm = RoiManager.getRoiManager()
  rm.runCommand("Associate", "true")
  rm.runCommand("Show All without labels")
  #This fits ellipses to the nuclear mask
  IJ.run(nucMaskIp, "Ellipse Split", "binary=[Use standard watershed] add_to_manager add_to_results_table merge_when_relativ_overlap_larger_than_threshold overlap=95 major=0-Infinity minor=0-Infinity aspect=1-Infinity stack")
  nucMaskStack2 = ImageStack.create(nucIp.width, nucIp.height, nucStack.size(), cmStack.getProcessor(1).getBitDepth())
  rois = rm.getRoisAsArray()
  #This loop is used to setup an image for the voronoi diagram nucelar splitting, setting pixels which
  #are contested between two neighboring nuclei to 0
  for roi in rois:
    z = roi.getZPosition()- 1
    try:
	    for p in roi.getContainedPoints():
	      try:
	        val = nucMaskStack2.getVoxel(p.x, p.y, z)
	        if val == 1: pass
	        elif val == 2:
	          nucMaskStack2.setVoxel(p.x, p.y, z, 1)
	        else: nucMaskStack2.setVoxel(p.x, p.y, z, 2)
	      except: pass
    except: pass
  rm.reset()
  for proc in (nucMaskStack2.getProcessor(n) for n in range(1, nucMaskStack2.size() + 1)): proc.threshold(1)
  nucMaskIp2 = ImagePlus("Masked Nuclei", nucMaskIp.getStack().getProcessor(1))
  nucMaskIp2.setStack(nucMaskStack2)
  #This calculates the voronoi diagram
  IJ.run(nucMaskIp2, "8-bit", "")
  nucMaskStack2 = nucMaskIp2.getStack()
  IJ.run(nucMaskIp2, "Close", "stack")
  IJ.run(nucMaskIp2, "Invert", "stack")
  IJ.run(nucMaskIp2, "Voronoi", "stack")
  for proc in (nucMaskStack2.getProcessor(n) for n in range(1, nucMaskStack2.size() + 1)): proc.threshold(0)
  IJ.run(nucMaskIp2, "Invert", "stack")
  #This burns the boundary lines of the voronoi diagram into the nuclear mask image, completing the segmentation
  nucAnalyzeIp = ImageCalculator().run("AND create stack", nucMaskIp, nucMaskIp2)
  nucAnalyzeIp.show()
  #This measures the nuclei
  IJ.run("Set Measurements...", "area mean centroid fit shape feret's stack decimal=3")
  Analyzer.setRedirectImage(nucIpForMeasure)
  IJ.run(nucAnalyzeIp, "Analyze Particles...", "display stack")
  rt = ResultsTable.getResultsTable()
  groupString = '_'.join("%s=%s" % (key,''.join(val)) for (key,val) in zip(groupBy, outerPairs))
  rt.saveAs(outputFolder +r"/nuclei_" + groupString + ".csv")
  IJ.run("Clear Results")
  rt = ResultsTable.getResultsTable()
  #This links nuclei to their respective cardiomyocyte (if present). The cardiomyocyte count mask image
  #has each pixel within a given cardiomyocyte set to its id number, so the minimum value for a nucleus roi in the cardiomyocyte
  #count mask image will be that cardiomyocyte's id number if the nucleus is 100% contained within a cardiomyocyte
  #or zero if not
  Analyzer.setRedirectImage(cmIp)
  IJ.run("Set Measurements...", "min decimal=3")
  IJ.run(nucAnalyzeIp, "Analyze Particles...", "display stack")
  rt.saveAs(outputFolder +r"/nucleilink_" + groupString + ".csv")
  IJ.run("Clear Results")
  rt = ResultsTable.getResultsTable()
  #This measures the cardiomyocytes
  Analyzer.setRedirectImage(cmIp)
  IJ.run("Set Measurements...", "area mean centroid shape feret's stack decimal=3")
  IJ.run(cmMaskIp, "Analyze Particles...", "size=500-Infinity display stack")
  rt.saveAs(outputFolder +r"/cm_" + groupString + ".csv")
  IJ.run("Clear Results")
  #This closes all open images, preventing the program from consuming more and more memory over time.
  IJ.run("Close All")
