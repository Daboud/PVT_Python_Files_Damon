# -*- coding: utf-8 -*-
"""
@author: damon.aboud
"""

import pvt_tools
import pvt_tools_5D
import matplotlib.pyplot as plt
import math

#--------------------------PROCESSING PARAMETERS-----------------------------
xmove=0.
scanv=1.0
xlength0=0.4
ylength0=0.4
holewidth=0.100
spacing=0.100
depthcount=1
Angle=math.pi/4
fatchange=.035
tallchange=0.020
LW=0.0109
ol0=0.0
zdrill=0.030
set_stage = '5D'
axo=1 #Axis orientation: Set 1 if you want to draw lines along +X and then raster-shift in the -Y direction (like reading a page of a book, right then down). 
#                        Set 2 if you want to draw lines along +Y and then raster-shift in the +X direction (up then right)
#-----------------------------------------------------------------------------


sa=math.sin(Angle)
ca=math.cos(Angle)
LWa=LW/ca
maxaccel=50.0
Zmaxaccel=0.8
zdepth=zdrill*(depthcount-1)
shuttertime=0.001

pitch=holewidth+spacing
zpitch=pitch*sa
ypitch=pitch*ca #not used 
Rtot=holewidth-LWa+tallchange
Rshift0=LWa*(1-ol0)
N0=Rtot/Rshift0
N=round(N0)
Rshift=Rtot/N
negRshifty=Rshift*ca * -1
Rshiftz=Rshift*sa
OLtrue=1-Rshift/LWa
# print(OLtrue)
N=int(N)

stripecount=round(xlength0/pitch)
stripecount=int(stripecount)
solidcount=round(ylength0/pitch)
solidcount=int(solidcount)

xlength=stripecount*pitch
ylength=solidcount*pitch

Droptime=math.sqrt(6*Rtot/maxaccel)
Elbowtime=4.*scanv/maxaccel
OSdroptime=math.sqrt(6*Rshiftz*N/Zmaxaccel)  #should only be used for full z change not single shift
Nextlinetime=math.sqrt(6*zpitch/Zmaxaccel)
shutterspace=0.001*scanv
zflytime=math.sqrt(6*Rshiftz/Zmaxaccel)
xflytime=3*(scanv**2+math.sqrt(scanv**2-(2/3*maxaccel)*-xlength)) / maxaccel

def format_pvt(pvt, set_stage, axo):
    if set_stage == '3D':
        return pvt
    elif set_stage == '5D':
        if axo==1:
            return pvt[0:7] + [0.0, 0.0, 0.0, 0.0] + pvt[7:9]
        if axo==2:
            return [pvt[0] , pvt[3]*-1 , pvt[4]*-1 , pvt[1] , pvt[2] , pvt[5] , pvt[6] , 0 , 0 , 0 , 0 , pvt[7] , pvt[8]]
    else:
        raise ('Error')
        
def stage_start_move():
    time=4*scanv/maxaccel #time calculated manually
    shutter=0
    shutter_sp=0
    x_p=0
    x_v=scanv
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def shutter_open(v):
    time=shuttertime #sec
    shutter=5
    shutter_sp=0
    x_v=v
    x_p=time*v
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def shutter_close(v):
    time=shuttertime #sec
    shutter=-5
    shutter_sp=-0
    x_v=v
    x_p=time*v
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def lase (x,x_vel):
    shutter=0
    shutter_sp=0
    x_p=x
    x_v=x_vel
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    time=(x_p/x_v) #sec
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def flyback (yp,zp):
    shutter=0
    shutter_sp=0
    x_p=-xlength
    x_v=scanv
    y_p=yp
    y_v=0
    z_p=zp
    z_v=0
    time=max(xflytime,zflytime)
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def left_dropdown (dy,dz):
    time=max(Elbowtime,Droptime,OSdroptime)
    shutter=0
    shutter_sp=0   #velocity tend to zero for both the z and y stages:
    x_p=0
    x_v=scanv
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def get_to_next_line ():
    shutter=0
    shutter_sp=0
    x_p=0
    x_v=0
    y_p=pitch*ca *-1
    y_v=0
    z_p=-1*zdrill*depthcount+pitch*sa
    z_v=0
    time=Nextlinetime #sec
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def change_places (xp,yp,zp):
    shutter=0
    shutter_sp=0
    x_p=xp
    x_v=0
    y_p=yp
    y_v=0
    z_p=zp
    z_v=0
    time=.001+max(math.sqrt(8*abs(xp)/maxaccel) , math.sqrt(8*abs(yp)/maxaccel) , math.sqrt(8*abs(zp)/Zmaxaccel))
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def initialize ():
    return format_pvt([1.23456789,0,0,0,0,0,0,0,0], set_stage, axo)
def write_PVT(PVT):
    #write a PVT file.txt
    f=open('da_SqHole_W'+str(round(holewidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt','w')
    for i in range (0,len(PVT)):
        f.write(" ".join(map(str, PVT[i])))
        f.write('\n')
    f.close()
def make_trajectory():

    PVT=[]
    PVT.append(initialize ())
    for l in range (solidcount):
        PVT.append(stage_start_move())
        for k in range (depthcount):
            for j in range (N):
                for i in range(stripecount): #lase several lines straight
                    PVT.append(shutter_open(scanv))
                    PVT.append(lase (holewidth-fatchange,scanv))
                    PVT.append(shutter_close(scanv))
                    PVT.append(lase (spacing+fatchange-shutterspace*2,scanv))
                PVT.append(flyback(negRshifty,Rshiftz))
            PVT.append(left_dropdown (-negRshifty*N, zdrill-N*Rshiftz))
        PVT.append(get_to_next_line ())
        
    PVT.append(change_places (0,pitch*ca*solidcount,-pitch*sa*solidcount)) #return to 0,0,0
    PVT.append(change_places (xmove,0,0)) #execute xmove

    write_PVT(PVT)
   
    xsum=sum(row[1] for row in PVT)
    ysum=sum(row[3] for row in PVT)
    zsum=sum(row[5] for row in PVT)
    print('Overlap % = ',str(format(OLtrue, '.4f')))
    print('Time (min) = ',str(format(sum(row[0] for row in PVT)/60, '.2f')))
    print('Sum X = ',str(format(xsum, '.5f')))
    print('Sum Y = ',str(format(ysum, '.5f')))
    print('Sum Z = ',str(format(zsum, '.5f')))
   
make_trajectory()


if set_stage == '3D':
    p = pvt_tools.PVT('da_SqHole_W'+str(round(holewidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt')
if set_stage == '5D':
    p = pvt_tools_5D.PVT('da_SqHole_W'+str(round(holewidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt')



p.read()
p.plot()


#plt.plot(p.A[:,1])
plt.plot(p.A[:,0])
plt.plot(p.A[:,1])
plt.plot(p.A[:,2])

plt.figure()
plt.plot(p.V[:,1]) #Y
plt.plot(p.V[:,2]) #Z
