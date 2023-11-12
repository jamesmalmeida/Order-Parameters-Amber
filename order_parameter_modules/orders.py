import numpy as np
from skspatial.objects import Plane
from skspatial.objects import Points

def watorderparam(l,molnum,molrange,coord,box,atm1,atm2):
	hc = np.zeros((3))
	dip = np.zeros((3))
	h1 = molrange[molnum][0]+atm1
	h2 = molrange[molnum][0]+atm2
	oid = molrange[molnum][0]
	hc[:] = (coord[l,h1,:]+coord[l,h2,:])/2 
	for k in range(0,3):
		hc[k] = (coord[l,h1,k]+coord[l,h2,k])/2
		if (coord[l,h1,k]-coord[l,h2,k]) > (box[l,k]/2):
			hc[k] = hc[k] - (box[l,k]/2)
		dip[k] = hc[k] - coord[l,oid,k]
		if dip[k] > (box[l,k]/2):
			dip[k] = dip[k] - box[l,k]/2
	dotp = np.dot(dip/np.linalg.norm(dip),[0,0,1])
	angle = np.arccos(dotp)
	orderparam = (0.5*(3*(np.cos(angle))**2-1))
	return orderparam

def orderparamring(l,molnum,molrange,coord,atm1,atm2,atm3):
	p1idx = molrange[molnum][0]+atm1
	p2idx = molrange[molnum][0]+atm2
	p3idx = molrange[molnum][0]+atm3
	p1 = coord[l,p1idx,:]
	p2 = coord[l,p2idx,:]
	p3 = coord[l,p3idx,:]
	v1 = p3 - p1
	v2 = p2 - p1
	normal = np.cross(v1, v2)
	dotp = np.dot(normal/np.linalg.norm(normal),[0,0,1])
	angle = np.arccos(dotp) #*57.2958
	order = (0.5*(3*(np.cos(angle))**2-1)) 
	return order

def calccurtosis(l,molnum,molrange,coord,atm1,atm2,refdist):
	p1idx = molrange[molnum][0]
	p2idx = molrange[molnum][0]+5
	p1 = coord[l,p1idx,:]
	p2 = coord[l,p2idx,:]
	dist = np.linalg.norm(p1-p2)
	curtosis = dist - refdist
	return curtosis

def orderparamlinear(l,molnum,molrange,coord,atm1,atm2):
	p1idx = molrange[molnum][0]+atm1
	p2idx = molrange[molnum][0]+atm2
	p1 = coord[l,p1idx,:]
	p2 = coord[l,p2idx,:]
	vec = p2-p1
	dotp = np.dot(vec/np.linalg.norm(vec),[0,0,1])
	angle = np.arccos(dotp)
	order = (0.5*(3*(np.cos(angle))**2-1))
	return order

def orderparamplane(l,molnum,molrange,coord,atm1,atm2):
	points = Points(coord[l,molrange[molnum][0]+47:molrange[molnum][0]+190:1,:])
	#print(points)
	#print(type(points))
	plane = Plane.best_fit(points)
	normal = plane.normal
	dotp = np.dot(normal/np.linalg.norm(normal),[0,0,1])
	angle = np.arccos(dotp)
	order = (0.5*(3*(np.cos(angle))**2-1))
	return order
