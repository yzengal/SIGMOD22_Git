import os
import sys
import random
import commands
import multiprocessing

class CFR:
	caseN = 10
	execNames = ["FRT", "Bar96", "Bar98", "DCnn", "DCsam"]

def run(execName, dim, pointFile, inputFile, logFile):
	cmdLine = "./%s-%d %s %s" % (execName, dim, pointFile, inputFile)
	print(cmdLine)
	line = commands.getoutput(cmdLine)
	print "[Finish]", cmdLine, logFile
	with open(logFile, "w") as fout:
		fout.write(line)


def batchRun(srcFilePath, desFilePath, nprocess):
	bid, eid = 0, CFR.caseN
	pool = multiprocessing.Pool(processes = nprocess)
	
	if not os.path.exists(desFilePath):
		os.mkdir(desFilePath)
	pointFiles = ["nor", "uni", "exp", "skew"]
	dimList = [2, 3, 4, 5, 10, 20, 100]
	alpha = 2
	for pointFilePrefix in pointFiles:
		setFilePath = os.path.join(desFilePath, pointFilePrefix)
		if not os.path.exists(setFilePath):
			os.mkdir(setFilePath)
		for dim in dimList:
			pointFileName = "%d.txt" % (dim)
			pointFile = os.path.join(srcFilePath, pointFilePrefix, pointFileName)
			for execName in CFR.execNames:
				execFilePath = os.path.join(setFilePath, execName)
				if not os.path.exists(execFilePath):
					os.mkdir(execFilePath)
				tmpFilePath = os.path.join(execFilePath, str(dim))
				if not os.path.exists(tmpFilePath):
					os.mkdir(tmpFilePath)				
				for caseId in xrange(bid, eid):
					inputFileName = "data_%02d.txt" % (caseId)
					inputFile = os.path.join(srcFilePath, str(alpha), inputFileName)
					desFileName = "data_%02d.txt" % (caseId)
					desFile = os.path.join(execFilePath, str(dim), desFileName)
					if os.path.exists(desFile):
						continue
					pool.apply_async(run, (execName, dim, pointFile, inputFile, desFile, ))

	pool.close()
	pool.join()

	
def exp():
	nprocess = 1
	srcFilePath = "./multiData"
	desFilePath = "./multiResult"
	batchRun(srcFilePath, desFilePath, nprocess)

	
if __name__ == '__main__':
	exp()
