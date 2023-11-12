from parmed.amber import AmberParm
from parmed.amber import netcdffiles
import numpy as np
import time
import sys
from order_parameter_modules.readpdb import readpdbfile
from order_parameter_modules.orders import watorderparam, orderparamring, calccurtosis, orderparamlinear, orderparamplane 

#inputtrajectory = 'npt_prod2.nc'
inputtrajectory = 'npt_prod.nc'
#inputtrajectory = 'nvt_init.nc'

#READING NETCDF FILE
parm = AmberParm('system.prmtop',xyz='system.inpcrd')
mass = (np.array(parm.parm_data["MASS"]))/6.02E23
charge = (np.array(parm.parm_data["CHARGE"]))
natoms = len(parm.atoms)
print('Found',natoms,'atoms')
coord = netcdffiles.NetCDFTraj.open_old(inputtrajectory).coordinates
nstep = np.ma.size(coord,0)
print('Found',nstep,'Steps')
box = netcdffiles.NetCDFTraj.open_old(inputtrajectory).box
boxz = box[0,2]
bins = int(boxz)
print('Found',bins,'bins')

#READING PDB
residue, atomtype, groupname, resnum, uniqueres = readpdbfile('system.pdb')

#MANUALLY SETTING THE NUMBER OF ATOMS PER MOLECULE
atmmolnumber = {
	"CTB": 62,
	"CL": 1,
	"POL": 191,
	"HEX": 20,
	"HEP": 23,
	"OCT": 26,
	"NON": 29,
	"CYX": 18,
	"CYP": 21,
	"BEN": 12,
	"TOL": 15,
	"WAT": 3
}

#IMPROVE THIS ROUTINE TO AUTOMATICALLY OBTAIN THE NUMBER OF ATOMS PER MOLECULE RESIDUE
#moleculesperres = {}
#for res in uniqueres:
#    atmnumres = int(residue.count(res))
#    molnumres = int(residue.count(res)/atmmolnumber[res])
#    moleculesperres.update({res: molnumres})
#
#print(moleculesperres)
#input("Press Enter to continue...")

fullorderparambin = []
fullanglebenzenebin = []
fullangletoluenebin = []
fullanglecyclohexanebin = []
fullanglecycloheptanebin = []
fullcurtosishexbin = []
fullcurtosishepbin = [] 
fullcurtosisoctbin = []
fullcurtosisnonbin = []
fullcurtosisctbbin = []
fullanglectbbin = []
fullcurtosispolbin = []
fullanglepolbin = []

for l in range(nstep):
    start_time = time.time()
    #CALCULATION THE CENTER OF MASS
    cm = []
    molnum = -1
    molrange = []
    molres = []
    i = 0
    while i <= natoms: 
        start = i
        try:
            end = i + atmmolnumber[residue[i]] - 1
        except:
            break
        molnum += 1
        cm.append(0)
        molmass = 0
        molrange.append([start, end])
        molres.append(residue[i])
        i = i + atmmolnumber[residue[i]]
        for j in range(start,end+1):
            cm[molnum] = cm[molnum] + coord[l,j,2]*mass[j]
            molmass = molmass + mass[j]
        cm[molnum] = cm[molnum]/molmass
    
    totalmolnumber = molnum #CHANGE ALL THE OTHER totalmolnumber to molnum, to get rid of this variable
   
    #Initializing all the np arrays for the order parameter calculations
    orderparam = np.zeros((totalmolnumber))
    anglebenzene = np.zeros(totalmolnumber)
    angletoluene = np.zeros(totalmolnumber)
    anglecyclohexane = np.zeros(totalmolnumber)
    anglecycloheptane = np.zeros(totalmolnumber)
    curtosishex = np.zeros(totalmolnumber)
    curtosishep = np.zeros(totalmolnumber)
    curtosisoct = np.zeros(totalmolnumber)
    curtosisnon = np.zeros(totalmolnumber)
    curtosisctb = np.zeros(totalmolnumber)
    anglectb = np.zeros(totalmolnumber)
    curtosispol = np.zeros(totalmolnumber)
    anglepol = np.zeros(totalmolnumber)
    #Calculate the order parameters, could paralelize this loop
    for molnum in range(0,totalmolnumber):
        if molres[molnum] == 'WAT':
            orderparam[molnum] = watorderparam(l,molnum,molrange,coord,box,1,2)
        if molres[molnum] == 'BEN':
            anglebenzene[molnum] = orderparamring(l,molnum,molrange,coord,0,2,4)
        if molres[molnum] == 'TOL':
            angletoluene[molnum] = orderparamring(l,molnum,molrange,coord,1,3,5)
        if molres[molnum] == 'CYX':
            anglecyclohexane[molnum] = orderparamring(l,molnum,molrange,coord,0,2,4)
        if molres[molnum] == 'CYP':
            anglecycloheptane[molnum] = orderparamring(l,molnum,molrange,coord,0,2,5)
        if molres[molnum] == 'HEX':
            curtosishex[molnum] = calccurtosis(l,molnum,molrange,coord,0,5,6.435408)
        if molres[molnum] == 'HEP':
            curtosishep[molnum] = calccurtosis(l,molnum,molrange,coord,0,6,7.656316)
        if molres[molnum] == 'OCT':
            curtosisoct[molnum] = calccurtosis(l,molnum,molrange,coord,0,7,8.975528)
        if molres[molnum] == 'NON':
            curtosisoct[molnum] = calccurtosis(l,molnum,molrange,coord,0,8,10.213422)
        if molres[molnum] == 'CTB':
            curtosisctb[molnum] = calccurtosis(l,molnum,molrange,coord,0,58,20.450132)
            anglectb[molnum] = orderparamlinear(l,molnum,molrange,coord,0,4)
        if molres[molnum] == 'POL':
            curtosispol[molnum] = calccurtosis(l,molnum,molrange,coord,2,40,16.386875)
            anglepol[molnum] = orderparamplane(l,molnum,molrange,coord,47,190)
    #INITIALIZING THE BIN RESOLVED ORDER PARAMETER VARIABLES
    orderparambin = np.zeros((bins))
    anglebenzenebin = np.zeros((bins))
    angletoluenebin = np.zeros((bins))
    anglecyclohexanebin = np.zeros((bins))
    anglecycloheptanebin = np.zeros((bins))
    curtosishexbin = np.zeros((bins))
    curtosishepbin = np.zeros((bins))
    curtosisoctbin = np.zeros((bins))
    curtosisnonbin = np.zeros((bins))
    curtosisctbbin = np.zeros((bins))
    anglectbbin = np.zeros((bins))
    curtosispolbin = np.zeros((bins))
    anglepolbin = np.zeros((bins))
    #CALCULATE THE BIN RESOLVED ORDER PARAMETERS
    for k in range(0,bins):
        #URGENTLY NEEDS TO IMPLEMENT THE CASE SELECT FROM PYTHON 3.10
        #WATER
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'WAT'))[0] 
        orderparamacc = 0
        if len(rows) != 0:
            for j in rows:
                orderparamacc = orderparamacc + orderparam[j]
            orderparambin[k] = orderparamacc/len(rows)
        #BENZENE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'BEN'))[0]
        anglebenzeneacc = 0
        if len(rows) != 0:
            for j in rows:
                anglebenzeneacc = anglebenzeneacc + anglebenzene[j]
            anglebenzenebin[k] = anglebenzeneacc/len(rows)
        #TOLUENE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'TOL'))[0]
        angletolueneacc = 0
        if len(rows) != 0:
            for j in rows:
                angletolueneacc = angletolueneacc + angletoluene[j]
            angletoluenebin[k] = angletolueneacc/len(rows)
        #CYCLOHEXANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'CYX'))[0]
        anglecyclohexaneacc = 0
        if len(rows) != 0:
            for j in rows:
                anglecyclohexaneacc = anglecyclohexaneacc + anglecyclohexane[j]
            anglecyclohexanebin[k] = anglecyclohexaneacc/len(rows)
        #CYCLOHEPTANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'CYP'))[0]
        anglecycloheptaneacc = 0
        if len(rows) != 0:
            for j in rows:
                anglecycloheptaneacc = anglecycloheptaneacc + anglecycloheptane[j]
            anglecycloheptanebin[k] = anglecycloheptaneacc/len(rows)
        #HEXANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'HEX'))[0]
        curtosishexacc = 0
        if len(rows) != 0:
            for j in rows:
                curtosishexacc = curtosishexacc + curtosishex[j]
            curtosishexbin[k] = curtosishexacc/len(rows)
        #HEPTANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'HEP'))[0]
        curtosishepacc = 0
        if len(rows) != 0:
            for j in rows:
                curtosishepacc = curtosishepacc + curtosishep[j]
            curtosishepbin[k] = curtosishepacc/len(rows)
        #OCTANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'OCT'))[0]
        curtosisoctacc = 0
        if len(rows) != 0:
            for j in rows:
                curtosisoctacc = curtosisoctacc + curtosisoct[j]
            curtosisoctbin[k] = curtosisoctacc/len(rows)
        #NONANE
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'NON'))[0]
        curtosisnonacc = 0
        if len(rows) != 0:
            for j in rows:
                curtosisnonacc = curtosisnonacc + curtosisnon[j]
            curtosisnonbin[k] = curtosisnonacc/len(rows)
        #CTAB
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'CTB'))[0]
        curtosisctbacc = 0
        anglectbacc =0
        if len(rows) != 0:
            for j in rows:
                curtosisctbacc = curtosisctbacc + curtosisctb[j]
                anglectbacc = anglectbacc + anglectb[j]
            curtosisctbbin[k] = curtosisctbacc/len(rows)
            anglectbbin[k] = anglectbacc/len(rows)
        #POL
        rows = np.where((np.asarray(cm[0:totalmolnumber]) >= k) & (np.asarray(cm[0:totalmolnumber]) < k + 1) & (np.asarray(molres[0:totalmolnumber]) == 'POL'))[0]
        curtosispolacc = 0
        anglepolacc = 0
        if len(rows) != 0:
            for j in rows:
                curtosispolacc = curtosispolacc + curtosispol[j]
                anglepolacc = anglepolacc + anglepol[j]
            curtosispolbin[k] = curtosispolacc/len(rows)
            anglepolbin[k] = anglepolacc/len(rows)

    #APPEND TO TOTALS, URGENTLY NEED TO CHANGE THIS TO AN ARRAY
    fullorderparambin.append(orderparambin)
    fullanglebenzenebin.append(anglebenzenebin)
    fullangletoluenebin.append(angletoluenebin)
    fullanglecyclohexanebin.append(anglecyclohexanebin)
    fullanglecycloheptanebin.append(anglecycloheptanebin)
    fullcurtosishexbin.append(curtosishexbin)
    fullcurtosishepbin.append(curtosishepbin)
    fullcurtosisoctbin.append(curtosisoctbin)
    fullcurtosisnonbin.append(curtosisnonbin)
    fullcurtosisctbbin.append(curtosisctbbin)
    fullcurtosispolbin.append(curtosispolbin)
    fullanglectbbin.append(anglectbbin)
    fullanglepolbin.append(anglepolbin)
    print("Step",l,"time",(time.time() - start_time),'s')

#CALCULATE THE AVERAGE FOR EACH BIN
fullorderparambinavg = np.zeros((bins))
fullanglebenzenebinavg = np.zeros((bins))
fullangletoluenebinavg = np.zeros((bins))
fullanglecyclohexanebinavg = np.zeros((bins))
fullanglecycloheptanebinavg = np.zeros((bins))
fullcurtosishexbinavg = np.zeros(bins)
fullcurtosishepbinavg = np.zeros(bins)
fullcurtosisoctbinavg = np.zeros(bins)
fullcurtosisnonbinavg = np.zeros(bins)
fullcurtosisctbbinavg = np.zeros(bins)
fullcurtosispolbinavg = np.zeros(bins)
fullanglectbbinavg = np.zeros(bins)
fullanglepolbinavg = np.zeros(bins)

#print(len(fullorderparambin))
#print(len(fullorderparambin[0]))
#fullorderparambinarr = np.array(fullorderparambin)
#print(fullorderparambinarr.shape)

fullorderparambinavg        = np.average(np.array(fullorderparambin), axis=0)
fullanglebenzenebinavg      = np.average(np.array(fullanglebenzenebin), axis=0)
fullangletoluenebinavg      = np.average(np.array(fullangletoluenebin), axis=0)
fullanglecyclohexanebinavg  = np.average(np.array(fullanglecyclohexanebin), axis=0)
fullanglecycloheptanebinavg = np.average(np.array(fullanglecycloheptanebin), axis=0)
fullcurtosishexbinavg       = np.average(np.array(fullcurtosishexbin), axis=0)
fullcurtosishepbinavg       = np.average(np.array(fullcurtosishepbin), axis=0)
fullcurtosisoctbinavg       = np.average(np.array(fullcurtosisoctbin), axis=0)
fullcurtosisnonbinavg       = np.average(np.array(fullcurtosisnonbin), axis=0)
fullcurtosisctbbinavg       = np.average(np.array(fullcurtosisctbbin), axis=0)
fullcurtosispolbinavg       = np.average(np.array(fullcurtosispolbin), axis=0)
fullanglectbbinavg          = np.average(np.array(fullanglectbbin), axis=0)
fullanglepolbinavg          = np.average(np.array(fullanglepolbin), axis=0)
#OLD FORTRANISH WAY TO AVERAGE AND NORMALIZE
#for l in range(nstep):
#    for i in range(bins):
#        fullorderparambinavg[i] = fullorderparambinavg[i] + fullorderparambin[l][i]
#        fullanglebenzenebinavg[i] = fullanglebenzenebinavg[i] + fullanglebenzenebin[l][i]
#        fullangletoluenebinavg[i] = fullangletoluenebinavg[i] + fullangletoluenebin[l][i]
#        fullanglecyclohexanebinavg[i] = fullanglecyclohexanebinavg[i] + fullanglecyclohexanebin[l][i]
#        fullanglecycloheptanebinavg[i] = fullanglecycloheptanebinavg[i] + fullanglecycloheptanebin[l][i]
#        fullcurtosishexbinavg[i] = fullcurtosishexbinavg[i] + fullcurtosishexbin[l][i]
#        fullcurtosishepbinavg[i] = fullcurtosishepbinavg[i] + fullcurtosishepbin[l][i]
#        fullcurtosisoctbinavg[i] = fullcurtosisoctbinavg[i] + fullcurtosisoctbin[l][i]
#        fullcurtosisnonbinavg[i] = fullcurtosisnonbinavg[i] + fullcurtosisnonbin[l][i]
#        fullcurtosisctbbinavg[i] = fullcurtosisctbbinavg[i] + fullcurtosisctbbin[l][i]
#        fullcurtosispolbinavg[i] = fullcurtosispolbinavg[i] + fullcurtosispolbin[l][i]
#        fullanglectbbinavg[i] = fullanglectbbinavg[i] + fullanglectbbin[l][i]
#        fullanglepolbinavg[i] = fullanglepolbinavg[i] + fullanglepolbin[l][i]
#fullorderparambinavg = fullorderparambinavg/nstep
#fullanglebenzenebinavg = fullanglebenzenebinavg/nstep
#fullangletoluenebinavg = fullangletoluenebinavg/nstep
#fullanglecyclohexanebinavg = fullanglecyclohexanebinavg/nstep
#fullanglecycloheptanebinavg = fullanglecycloheptanebinavg/nstep
#fullcurtosishexbinavg = fullcurtosishexbinavg/nstep
#fullcurtosishepbinavg = fullcurtosishepbinavg/nstep
#fullcurtosisoctbinavg = fullcurtosisoctbinavg/nstep
#fullcurtosisnonbinavg = fullcurtosisnonbinavg/nstep
#fullcurtosisctbbinavg = fullcurtosisctbbinavg/nstep
#fullcurtosispolbinavg = fullcurtosispolbinavg/nstep
#fullanglectbbinavg = fullanglectbbinavg/nstep
#fullanglepolbinavg = fullanglepolbinavg/nstep

x = np.linspace(0,bins,bins)
lenght=len(x)
z = np.column_stack((x,np.asarray(fullorderparambinavg)))
np.savetxt("order-parameter.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullanglebenzenebinavg)))
np.savetxt("order-parameter-benzene.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullangletoluenebinavg)))
np.savetxt("order-parameter-toluene.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullanglecyclohexanebinavg)))
np.savetxt("order-parameter-cyclohexane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullanglecycloheptanebinavg)))
np.savetxt("order-parameter-cycloheptane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosishexbinavg)))
np.savetxt("order-parameter-hexane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosishepbinavg)))
np.savetxt("order-parameter-heptane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosisoctbinavg)))
np.savetxt("order-parameter-octane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosisnonbinavg)))
np.savetxt("order-parameter-nonane.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosisctbbinavg)))
np.savetxt("order-parameter-ctab.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullcurtosispolbinavg)))
np.savetxt("order-parameter-hpg11.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullanglectbbinavg)))
np.savetxt("order-parameter-angle-ctab.dat", z, delimiter=" ")
z = np.column_stack((x,np.asarray(fullanglepolbinavg)))
np.savetxt("order-parameter-angle-hpg11.dat", z, delimiter=" ")

