import os
import sys
from random import randint, sample, shuffle
from math import log


class constForReal():
	caseN = 10
	pointFiles = ["checkinNYC.txt", "checkinTKY.txt", "gaiaChengdu.txt", "gaiaHaikou.txt"]
	alphas = [2]
	
class CFR(constForReal):
	pass
	
def genBeta(alpha, n):
	lnk = log(alpha, 10)
	base, beg = 1000, 1
	dx = (alpha - beg) * 1.0 / 1000;
	P = [0.0] * base;
	sum = 0.0;
	
	for i in xrange(base):
		P[i] = dx / (lnk * (beg + i*dx))
		sum += P[i]
	
	for i in xrange(base):
		P[i] /= sum;
		if i>0: 
			P[i] += P[i-1]
		
	
	retList = []
	for j in xrange(n):
		ret = alpha;
		randf = randint(0, base*10-1) / (base*10.0)
		for i in xrange(base):
			if randf <= P[i]:
				ret = beg + i * dx
				break
		ret /= alpha
		ret = min(ret, 1.0)
		ret = max(1.0/alpha, ret)
		retList.append(ret)
		
	return retList
	
def genRealData(desPath, caseN):
	for k in CFR.alphas:
		tmpFilePath = os.path.join(desPath, str(k))
		if not os.path.exists(tmpFilePath):
			os.mkdir(tmpFilePath)		
		betas = genBeta(k, caseN)
		for i in xrange(caseN):
			desFileName = "data_%02d.txt" % (i)
			desFileName = os.path.join(tmpFilePath, desFileName)
			with open(desFileName, "w") as fout:
				fout.write("%.4f %d\n" % (betas[i], k))
	
	
def genRealDataSet(desPath, caseN):
	if not os.path.exists(desPath):
		os.mkdir(desPath)
	genRealData(desPath, caseN)
	
	
def exp0():
	desPath = "realData"
	genRealDataSet(desPath, CFR.caseN)
	
	
if __name__ == "__main__":
	exp0()
	