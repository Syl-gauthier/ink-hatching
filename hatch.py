import math
from types import *
from datetime import datetime
import sys
sys.path.append('/usr/share/inkscape/extensions')
import inkex

class Point(object):
	def __init__(self, XCoordinate, YCoordinate):
		if type(XCoordinate) is not FloatType:
			raise TypeError("Type of XCoordinate is not float")
		elif type(YCoordinate) is not FloatType:
			raise TypeError("Type of YCoordinate is not float")
		self.x = XCoordinate
		self.y = YCoordinate

	def __str__(self):
		return str(self.x) + "," + str(self.y)

	def __repr__(self):
		return str(self)

	def __add__(self, point):
		return Point(self.x + point.x, self.y + point.y)


class CurvePoint(Point):
	def __init__(self, TCoordinate, XCoordinate, YCoordinate):
		if type(TCoordinate) is not FloatType:
			raise TypeError("Type of TCoordinate is not float")
		elif type(XCoordinate) is not FloatType:
			raise TypeError("Type of XCoordinate is not float")
		elif type(YCoordinate) is not FloatType:
			raise TypeError("Type of YCoordinate is not float")
		self.t = TCoordinate
		self.x = XCoordinate
		self.y = YCoordinate


class Line(object):
	def __init__(self, slope, yIntercept=0.0):
		if type(slope) is not FloatType:
			raise TypeError("Type of slope should be float")
		elif type(yIntercept) is not FloatType:
			raise TypeError("Type of yIntercept should be float")
		self.a = slope
		self.b = yIntercept

	def __str__(self):
		firstPoint = Point(0.0, self.b)
		secondY = self.a*744.0 + self.b
		secondPoint = Point(744.0, secondY)
		return "M " + str(firstPoint) + " " + str(secondPoint)

	def __repr__(self):
		return "line{a:" + self.a + ",b:" + self.b + "}"


class Curve(object):
	def __init__(self, points):
		for point in points:
			if not isinstance(point, Point):
				raise TypeError("invalid point instance")
		self.curveType = "undefined Curve"
	
	def getSegments(self):
		raise AbstractFunctionError("This function is not callable")

	def getPoint(t):
		raise AbstractFunctionError("This function is not callable")

	def findIntersection(self, line):
		segments = self.getApprox()
		return segments.findIntersection(line)

	def quickIntersection(self, line):
		for segment in self.getSegments():
			if segment.quickIntersection(line):
				return True
		return False

	def getRoughLength(self):
		raise AbstractFunctionError("This function is not callable")

	def getApprox(self):
		length = self.getRoughLength()
		step = 1.0 / length
		
		segments = SegmentedCurve(self.getPoint(0.0), self.getPoint(1.0))
		for i in range(1, int(length)+1):
			segments.addSegment(self.getPoint(i*step))
		segments.addSegment(segments.lastPoint)
		
		return segments

	def split(self):
		raise AbstractFunctionError("This function is not callable")

	def extract(self, t1, t2):
		firstSplit = self.split(t2)[0]
		t1 = t1/t2
		secondSplit = firstSplit.split(t1)[1]
		return secondSplit


class SegmentCurve(Curve):
	def __init__(self, points):
		super(SegmentCurve, self).__init__(points)
		self.p0 = points[0]
		self.p1 = points[1]
		self.curveType = "segment"
	
	def __str__(self):
		return "M "+str(self.p0)+" "+str(self.p1)

	def getSegments(self):
		return [self]

	def getPoint(self, t):
		return CurvePoint(t, self.p0.x*(1-t)+self.p1.x*t, self.p0.y*(1-t)+self.p1.y*t)

	def findIntersection(self, line):
		denom = line.a * (self.p1.x - self.p0.x) - self.p1.y + self.p0.y
		if (denom==0):
			return False
		t = (self.p0.y - line.a * self.p0.x - line.b) / (line.a * (self.p1.x - self.p0.x) - self.p1.y + self.p0.y)


		if (not(0<=t<=1)):
			return False

		else:
			x = self.p0.x * (1 - t) + self.p1.x * t
			y = self.p0.y * (1 - t) + self.p1.y * t
			
			point = {"t": t, "x": x}
			return([point])

	def quickIntersection(self, line):
		denom = line.a * (self.p1.x - self.p0.x) - self.p1.y + self.p0.y
		if (denom==0):
			return False
		t = (self.p0.y - line.a * self.p0.x - line.b) / denom
		if (0<=t<=1) :
			return True

		return False

	def getRoughLength(self): #exact length
		return (math.sqrt((self.p1.x-self.p0.x)**2 + (self.p1.y-self.p0.y)**2))

	def split(self, t):
		curve1 = SegmentCurve([self.p0, self.getPoint(t)])
		curve2 = SegmentCurve([self.getPoint(t), self.p1])
		return [curve1, curve2]


class QuadraticCurve(Curve):
	def __init__(self, points):
		super(QuadraticCurve, self).__init__(points)
		self.p0 = points[0]
		self.p1 = points[1]
		self.p2 = points[2]

	def __str__(self):
		return "M "+str(self.p0)+" Q "+str(self.p1)+" "+str(self.p2)

	def getSegments(self):
		segment1 = SegmentCurve([self.p0, self.p1])
		segment2 = SegmentCurve([self.p1, self.p2])
		return [segment1, segment2]

	def getPoint(self, t):
		x = self.p0.x*(1-t)**2 + 2*self.p1.x*(1-t)*t + self.p2.x*t**2
		y = self.p0.y*(1-t)**2 + 2*self.p1.y*(1-t)*t + self.p2.y*t**2
		return CurvePoint(t, x, y)

	def getRoughLength(self):
		l1 = math.sqrt((self.p1.x-self.p0.x)**2 + (self.p1.y-self.p0.y)**2)
		l2 = math.sqrt((self.p2.x-self.p1.x)**2 + (self.p2.y-self.p1.y)**2)
		return l1+l2

	def split(self, t):
		pt = self.getPoint(t)
		pCurve1 = SegmentCurve([self.p0, self.p1]).getPoint(t)
		pCurve2 = SegmentCurve([self.p1, self.p2]).getPoint(t)
		curve1 = QuadraticCurve([self.p0, pCurve1, pt])
		curve2 = QuadraticCurve([pt, pCurve2, self.p2])
		return [curve1, curve2]


class CubicCurve(Curve):
	def __init__(self, points):
		super(CubicCurve, self).__init__(points)
		self.p0 = points[0]
		self.p1 = points[1]
		self.p2 = points[2]
		self.p3 = points[3]

	def __str__(self):
		return "M "+str(self.p0)+" C "+str(self.p1)+" "+str(self.p2)+" "+str(self.p3)

	def getSegments(self):
		segment1 = SegmentCurve([self.p0, self.p1])
		segment2 = SegmentCurve([self.p1, self.p2])
		segment3 = SegmentCurve([self.p2, self.p3])

		return [segment1, segment2, segment3]

	def getPoint(self, t):
		x = self.p0.x*(1-t)**3 + 3*self.p1.x*t*(1-t)**2 + 3*self.p2.x*(1-t)*t**2 + self.p3.x*t**3
		y = self.p0.y*(1-t)**3 + 3*self.p1.y*t*(1-t)**2 + 3*self.p2.y*(1-t)*t**2 + self.p3.y*t**3
		return CurvePoint(t, x, y)


	def getRoughLength(self):
		l1 = math.sqrt((self.p1.x-self.p0.x)**2 + (self.p1.y-self.p0.y)**2)
		l2 = math.sqrt((self.p2.x-self.p1.x)**2 + (self.p2.y-self.p1.y)**2)
		l3 = math.sqrt((self.p3.x-self.p2.x)**2 + (self.p3.y-self.p2.y)**2)
		return l1+l2+l3

	def split(self, t):
		p01 = SegmentCurve([self.p0, self.p1]).getPoint(t)
		p02 = QuadraticCurve([self.p0, self.p1, self.p2]).getPoint(t)
		p03 = self.getPoint(t)
		p12 = QuadraticCurve([self.p1, self.p2, self.p3]).getPoint(t)
		p21 = SegmentCurve([self.p2, self.p3]).getPoint(t)

		curve1 = CubicCurve([self.p0, p01, p02, p03])
		curve2 = CubicCurve([p03, p12, p21, self.p3])

		return [curve1, curve2]


class SegmentedCurve:
	def __init__(self, firstPoint, lastPoint):
		if not isinstance(firstPoint, CurvePoint) or not isinstance(lastPoint, CurvePoint):
			raise TypeError("SegmentedCurve firstPoint and lastPoint should be CurvePoints")
		self.firstPoint = firstPoint
		self.lastPoint = lastPoint
		self.segments = []

	def addSegment(self, curvePoint):
		if self.segments == []:
			newSegment = SegmentCurve([self.firstPoint, curvePoint]) 
		else:
			newSegment = SegmentCurve([self.segments[-1].p1, curvePoint])
		
		self.segments.extend([newSegment])
		return 0

	def findIntersection(self, line):
		intersections = []
		for segment in self.segments:
			intersection = segment.findIntersection(line)
			if intersection:
				t = segment.p0.t + intersection[0]["t"] * (segment.p1.t - segment.p0.t)
				intersections.append({"t": t, "x": intersection[0]["x"]})
		if intersections == []:
			return False
		else:
			return intersections


class Path(object):
	def __init__(self, SVGString):
		if type(SVGString) is not StringType:
			raise TypeError("SVGString should be a string")
		self.curves = []
		SVGString = SVGString.split()

		mode = "x"
		curvesList = []
		pathComplete = False
		curveIndex = 0
		for data in SVGString :
			if data in "mMcClLhHvVsSqQtTaAzZ" :
				mode = data

				if mode in "mM" :
					curveIndex = 1
					#print "path begin"
				elif (mode == "l") :
					#print "segment curve (rel)"
					curveIndex = 1
				elif (mode == "L") :
					#print "segment curve (abs)"
					curveIndex = 1
				elif (mode == "q") :
					#print "quadratic curve (rel)"
					curveIndex = 1
				elif (mode == "Q") :
					#print "quadratic curve (abs)"
					curveIndex = 1
				elif (mode == "c") :
					#print "cubic curve (rel)"
					curveIndex = 1
				elif (mode == "C") :
					#print "cubic curve (abs)"
					curveIndex = 1
				elif (mode in "zZ") :
					if len(curvesList[-1]) == 1:
						curvesList[-1].append(self.originPoint)
					pathComplete = True
				else :
					raise ValueError("undefined mode found in SVGString")

			else :
				coords = data.split(",")
				splitPoint = Point(float(coords[0]), float(coords[1]))
				if mode in "mM" :
					if curvesList == []:
						self.originPoint = splitPoint
					else:
						curvesList[-1].append(splitPoint)
					curvesList.append([splitPoint])
					
				elif mode == "l":
					basePoint = curvesList[-1][0]
					curvesList[-1].append(basePoint+splitPoint)
					curvesList.append([basePoint+splitPoint])
				elif mode == "L" :
					curvesList[-1].append(splitPoint)
					curvesList.append([splitPoint])
				elif mode == "q":
					if curveIndex%2 == 0 :
						basePoint = curvesList[-1][0]
						curvesList[-1].append(basePoint+splitPoint)
						curvesList.append([basePoint+splitPoint])
					else :
						basePoint = curvesList[-1][0]
						curvesList[-1].append(basePoint+splitPoint)
					curveIndex += 1
				elif mode == "Q" :
					if curveIndex%2 == 0 :
						curvesList[-1].append(splitPoint)
						curvesList.append([splitPoint])
					else :
						curvesList[-1].append(splitPoint)
					curveIndex += 1
				elif mode == "c" :
					if curveIndex%3 == 0 :
						basePoint = curvesList[-1][0]
						curvesList[-1].append(basePoint+splitPoint)
						curvesList.append([basePoint+splitPoint])
					else :
						basePoint = curvesList[-1][0]
						curvesList[-1].append(basePoint+splitPoint)
					curveIndex += 1
				elif mode == "C" :
					if curveIndex%3 == 0 :
						curvesList[-1].append(splitPoint)
						curvesList.append([splitPoint])
					else :
						curvesList[-1].append(splitPoint)
					curveIndex += 1

		for curve in curvesList:
			if len(curve) == 2:
				self.curves.extend([SegmentCurve(curve)])
			elif len(curve) == 3:
				self.curves.extend([QuadraticCurve(curve)])
			elif len(curve) == 4:
				self.curves.extend([CubicCurve(curve)])

	def setHatch(self, slope, width, offset=0.005):
		origin = self.originPoint.y - slope*self.originPoint.x + offset
		self.hatch = Hatch(slope, width, origin)

	def findIntersection(self, line):
		intersections = []
		for index, curve in enumerate(self.curves):
			intersection = curve.findIntersection(line)
			if intersection:
				for point in intersection:
					pathPoint = PathPoint(index, point["t"], point["x"])
					intersections.append(pathPoint)
		return intersections


class Hatch(object):
	def __init__(self, slope, width, origin):
		if type(slope) is not FloatType:
			raise TypeError("slope should be a float")
		if type(width) is not FloatType:
			raise TypeError("width should be a float")
		if type(origin) is not FloatType:
			raise TypeError("origin should be a float")
		self.slope = slope
		self.width = width
		self.origin = origin

	def getRange(self, curves):
		intersect = True
		maxIndex = 0
		minIndex = 0
		index = 0 

		firstIntersect = False
		#find a first valid line
		while firstIntersect == False:
			if index >=0:
				for curve in curves:
					if curve.quickIntersection(self.getLine(index)):
						firstIntersect = index
				if firstIntersect == False: 
					index += 1
				if index >= 10:
					index = -1
			elif index >= -10:
				for curve in curves:
					if curve.quickIntersection(self.getLine(index)):
						firstIntersect = index
				if firstIntersect == False: 
					index = index - 1
			else:
				inkex.debug("No intersection found")
				firstIntersect = True
		
		index = firstIntersect
		while(intersect == True):
			intersect = False
			for curve in curves:
				if curve.quickIntersection(self.getLine(index)):
					maxIndex = index+1
					intersect = True
			index += 1

		intersect = True
		index = firstIntersect
		while(intersect == True):
			intersect = False
			for curve in curves:
				if curve.quickIntersection(self.getLine(index)):
					minIndex = index -1
					intersect = True
			index -= 1
		return (minIndex, maxIndex)

	def getLine(self, index):
		return Line(self.slope, self.origin + index*self.width)

	def getHatchPath(self, path):

		def findPoint2(newPoints, point1):
			newPoints.sort(key = lambda point: point.x)
			for index, point in enumerate(newPoints):
				if point == point1:
					index1 = index
			if index1%2 == 0:
				return newPoints[index1+1]
			else:
				return newPoints[index1-1]

		def findClosest(curves, newPoints, segment, isLeft, hatchIndex):
			newPoints.sort()
			prevSegments = []
			for curve in curves:
				if curve[-1]["index"] == hatchIndex:
					prevSegments.append(curve[-1]["segment"])
			if segment[isLeft] > newPoints[-1]:
				nextPointIndex = 0
			else:
				for index, point in enumerate(newPoints):
					if point > segment[isLeft]:
						nextPointIndex = index
						break
			prevPointIndex = len(newPoints) -1 if nextPointIndex == 0 else nextPointIndex - 1
			if checkPoint([segment[isLeft], newPoints[nextPointIndex]], prevSegments):
				return (nextPointIndex, newPoints[nextPointIndex])
			elif checkPoint([newPoints[prevPointIndex], segment[isLeft]], prevSegments):
				return (prevPointIndex, newPoints[prevPointIndex])
			else:
				return (False, False)

		def checkPoint(limits, segments):
			for segment in segments:
				for point in segment:
					if isInLimit(limits, point):
						return False
			return True

		def isInLimit(limits, x):
			if x == limits[0] or x == limits[1]:
				return False
			elif limits[0] < limits[1]:
				if limits[0] < x < limits[1]:
					return True
				return False
			else:
				if limits[1] < x < limits[0]:
					return False
				return True	

		min, max = self.getRange(path.curves)
		curves = []
		for index in range(min, max):
			line = self.getLine(index)
			newPoints = path.findIntersection(line) #return PathPoint list
			for curve in curves:
				if curve[-1]["index"] == index-1:
					if newPoints != []:
						index1, point1 = findClosest(curves, newPoints, curve[-1]["segment"], 1, index-1)
						if point1 != False:
							point2 = findPoint2(newPoints, point1)
							newPoints.sort()
							for i, point in enumerate(newPoints):
								if point == point2:
									index2 = i
							if point2 != False:
								curve.append({"index": index, "segment": [point1, point2]})
								newPoints.sort()
								if index1 > index2:
									del newPoints[index1]
									del newPoints[index2]
								else:
									del newPoints[index2]
									del newPoints[index1]

			newPoints.sort(key= lambda point: point.x)
			for i in range(len(newPoints)/2):
				curves.append([{"index": index, "segment": [newPoints[2*i+1], newPoints[2*i]]}])
		return curves


class PathPoint:
	def __init__(self, index, t, x):
		if type(index) is not IntType:
			raise TypeError("PathPoint index should be an int")
		elif type(t) is not FloatType:
			raise TypeError("PathPoint t should be a Float")
		elif type(x) is not FloatType:
			raise TypeError("PathPoint x should be a Float")
		self.index = index
		self.t = t
		self.x = x
	
	def __str__(self):
		return "PathP{" + str(self.index) + ",t" + str(self.t) + "}"

	def __repr__(self):
		return str(self)

	def __eq__(self, point):
		if self.index == point.index and self.t == point.t:
			return True
		return False

	def __lt__(self, point):
		if self.index < point.index or self.index == point.index and self.t < point.t:
			return True
		return False

	def __le__(self, point):
		if self.index < point.index or self.index == point.index and self.t <= point.t:
			return True
		return False

	def __gt__(self, point):
		return not self <= point

	def __ge__(self, point):
		return not self < point


def getSVG(list, path):
	svgCurves = []

	for curve in list:
		string = ""
		for entry in curve:
			segment = entry["segment"]
			if string == "":
				string += "M "
				string += str(path.curves[segment[0].index].getPoint(segment[0].t))
				string += " L "
				string += str(path.curves[segment[1].index].getPoint(segment[1].t))
				string += " "
			else:
				string += str(path.curves[segment[0].index].getPoint(segment[0].t))
				string += " "
				string += str(path.curves[segment[1].index].getPoint(segment[1].t))
				string += " "
		svgCurves.append(string)
	return svgCurves


class HatchMaker(inkex.Effect):
	def __init__(self):
		inkex.Effect.__init__(self)
		
		self.OptionParser.add_option("", "--slope",
						action="store", type="float", 
						dest="slope", default="1.0",
						help="Slope of the hatches")
		self.OptionParser.add_option("", "--width",
						action="store", type="float", 
						dest="width", default="20.0",
						help="Width between two line")
		self.OptionParser.add_option("", "--offset",
						action="store", type="float", 
						dest="offset", default="0.0",
						help="Offset the hatches")

	def effect(self):		
		for id, node in self.selected.iteritems():
			if node.tag == inkex.addNS('path', 'svg'):
				d = node.get('d')				
				path = Path(d)
				radianAngle = self.options.slope*math.pi / 180
				slopeConverted = math.tan(radianAngle)
				widthConverted = self.options.width / math.cos(radianAngle)
				path.setHatch(slopeConverted, widthConverted, self.options.offset)
				list = path.hatch.getHatchPath(path)
				svgList = getSVG(list, path)
				for svg in svgList:
					inkex.etree.SubElement(self.current_layer, inkex.addNS("path", "svg"),{"d": svg, "style": "fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"})


'''
d ="m 101.01525,183.63101 l 289.91378,108.08633 288.90363,103.03556 -1.01015,-5.05076 -166.67517,226.27417 -171.72593,227.28429 -5.05077,1.0102 -97.9848,-90.91374 -97.9848,-90.91373 0,0 94.95434,-87.88327 89.90358,-97.9848 -119.198,-161.62441 -104.04571,-131.31984 0,0 z"
'''

if __name__ == '__main__':
	e = HatchMaker()
	e.affect()
