import os
import sys
from random import randint, sample, shuffle
import numpy as np
from math import log

class constForSyn():
	caseN = 10
	Range = 10**7
	mu = 0.5 * Range
	sigma = 0.10 * Range
	pointFiles = ["uni", "nor", "exp", "skew"]
	num = 5*10**4
	dimList = [2, 3, 4, 5, 10, 20, 100]
	alpha = 2
	
class CFR(constForSyn):
	pass

class baseGenerator:
	def gen(self, n):
		pass	
	
class randomGenerator(baseGenerator):

	def __init__(self, mx):
		self.mx = mx

	def setMx(self, mx):
		self.mx = mx

	def gen(self, n):
		ret = [0] * n
		for i in xrange(n):
			x = randint(-self.mx, self.mx)
			ret[i] = x
		return ret

class normalGenerator(baseGenerator):

	def __init__(self, mu, sigma):
		self.mu = mu
		self.sigma = sigma

	def gen(self, n, lb = None, rb = None):
		ret = np.random.normal(self.mu, self.sigma, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setMu(self, mu):
		self.mu = mu

	def setSigma(self, sigma):
		self.sigma = sigma


class uniformGenerator(baseGenerator):

	def __init__(self, low, high):
		self.low = low
		self.high = high

	def gen(self, n, lb = None, rb = None):
		ret = np.random.uniform(self.low, self.high, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setLow(self, low):
		self.low = low

	def setHigh(self, high):
		self.high = high

class expGenerator(baseGenerator):

	def __init__(self, mu):
		self.mu = mu

	def gen(self, n, lb = None, rb = None):
		ret = np.random.exponential(self.mu, n)
		for i in xrange(n):
			if lb is not None and ret[i]<lb:
				ret[i] = lb
			if rb is not None and ret[i]>rb:
				ret[i] = rb
		return ret

	def setMu(self, mu):
		self.mu = mu

def genPoints(gtor, n, dim):
	points = []
	vals = []
	for j in xrange(dim):
		vals.append(gtor.gen(n))
	mx = 0
	for j in xrange(dim):
		mx += (max(vals[j])-min(vals[j]))**2
	for i in xrange(n):
		tmp = []
		for j in xrange(dim):
			tmp.append(vals[j][i])
		points.append(tuple(tmp))
	ret = list(set(points))
	return ret
	
def genSkewedPoints(gtor, n, dim):
	alpha = 2
	points = []
	vals = []
	for j in xrange(dim):
		if j==0:
			vals.append(map(int, gtor.gen(n)))
		else:
			vals.append(map(lambda x:int(x)**alpha, gtor.gen(n)))
	mx = 0
	for j in xrange(dim):
		mx += (max(vals[j])-min(vals[j]))**2
	for i in xrange(n):
		tmp = []
		for j in xrange(dim):
			tmp.append(vals[j][i])
		points.append(tuple(tmp))
	ret = list(set(points))
	return ret
		

def genSynData(points, desPath, prefix, n, dim):
	tmpFilePath = os.path.join(desPath, prefix)
	if not os.path.exists(tmpFilePath):
		os.mkdir(tmpFilePath)
	desFileName = "%d.txt" % (dim)
	desFileName = os.path.join(tmpFilePath, desFileName)
	with open(desFileName, "w") as fout:
		fout.write("%d\n" % (len(points)))
		for i in xrange(len(points)):
			fout.write(" ".join(map(str, points[i]))+"\n")
	
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
	
def genParameter(desPath):
	k = CFR.alpha
	tmpFilePath = os.path.join(desPath, str(k))
	if not os.path.exists(tmpFilePath):
		os.mkdir(tmpFilePath)		
	betas = genBeta(k, CFR.caseN)
	for i in xrange(CFR.caseN):
		desFileName = "data_%02d.txt" % (i)
		desFileName = os.path.join(tmpFilePath, desFileName)
		with open(desFileName, "w") as fout:
			fout.write("%.4f %d\n" % (betas[i], k))	
		
def genSynDataSet(desPath, caseN):
	if not os.path.exists(desPath):
		os.mkdir(desPath)
	# gen data first
	for prefix in CFR.pointFiles:
		if prefix=="uni":
			mu = CFR.mu
			gtor = uniformGenerator(CFR.Range*(-1), CFR.Range)
			n = CFR.num
			for dim in CFR.dimList:
				points = genPoints(gtor, n, dim)
				genSynData(points, desPath, prefix, n, dim)
		elif prefix=="nor":
			mu, sig = CFR.mu, CFR.sigma
			gtor = normalGenerator(mu, sig)
			n = CFR.num
			for dim in CFR.dimList:
				points = genPoints(gtor, n, dim)
				genSynData(points, desPath, prefix, n, dim)		
		elif prefix=="exp":
			mu = CFR.mu
			gtor = expGenerator(mu)
			n = CFR.num
			for dim in CFR.dimList:
				points = genPoints(gtor, n, dim)
				genSynData(points, desPath, prefix, n, dim)	
		elif prefix=="skew":
			mu = CFR.mu
			gtor = uniformGenerator(CFR.Range*(-1), CFR.Range)
			n = CFR.num
			for dim in CFR.dimList:
				points = genSkewedPoints(gtor, n, dim)
				genSynData(points, desPath, prefix, n, dim)	
				
	# gen parameter beta second
	genParameter(desPath)
	
def exp0():
	desPath = "multiData"
	genSynDataSet(desPath, CFR.caseN)
	
	
if __name__ == "__main__":
	exp0()
	
