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
xlength0=0.5
ylength0=0.5
pwidth=0.100
spacing=0.070
depthcount=1
Angle=math.pi/4*44.5/45
fatchange= -0.025 #+ve makes pillar fatter
tallchange= -0.019#+ve makes pillars taller!
LW=0.0109
ol0=0.0
zdrill=0.030
set_stage = '5D'
axo=1 #Axis orientation: Set 1 if you want to draw lines along +X and then raster-shift in the -Y direction (like reading a page of a book, right then down). 
#                        Set 2 if you want to draw lines along +Y and then raster-shift in the +X direction (up then right)
#-----------------------------------------------------------------------------

shuttertime=0.001
sa=math.sin(Angle)
ca=math.cos(Angle)
LWa=LW/ca
maxaccel=50.0
Zmaxaccel=0.8
zdepth=zdrill*(depthcount-1)

pitch=pwidth+spacing
Rtot=pitch
Rshift0=LWa*(1-ol0)
N0=Rtot/Rshift0
N=round(N0)
Rshift=Rtot/N
Rshifty=Rshift*ca
Rshiftz=Rshift*sa
OLtrue=1-Rshift/LWa
N=int(N)

Ns=round((spacing-tallchange-LWa)/pitch*N)
Np=N-Ns
Ns=int(Ns)
Np=int(Np)

stripecount=round(xlength0/pitch)+1
solidcount=round(ylength0/pitch)
xlength=stripecount*pitch
ylength=solidcount*pitch


shutterspace=shuttertime*scanv
flydist=xlength-LW-pwidth-fatchange
flytime=(scanv+math.sqrt(scanv**2+4*(1/6*maxaccel)*(flydist))) / (2*1/6*maxaccel)

flydisto2=xlength/2
t2=( scanv + math.sqrt(scanv**2+2*maxaccel*flydisto2)) / maxaccel


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
    time=4*scanv/maxaccel #time calculated manually, a=8mm/s^2
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
    shutter_sp=0
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
def flyback (xp,yp,zp):
    shutter=0
    shutter_sp=0
    x_p=xp
    x_v=scanv
    y_p=yp
    y_v=0
    z_p=zp
    z_v=0
    time=max(flytime , math.sqrt(abs(6*yp/maxaccel)) , math.sqrt(abs(6*zp/Zmaxaccel)))
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def flyback1 (xp,yp,zp):
    shutter=0
    shutter_sp=0
    x_p=xp
    x_v=scanv-t2*maxaccel
    y_p=yp
    y_v=-Rshifty/t2
    z_p=zp
    z_v=Rshiftz/t2    #could cause errors
    time=t2
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def flyback2 (xp,yp,zp):
    shutter=0
    shutter_sp=0
    x_p=xp
    x_v=scanv
    y_p=yp
    y_v=0
    z_p=zp
    z_v=0
    time=t2
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def left_dropdown (dy,dz):
    time=max(4.*scanv/maxaccel , math.sqrt(6*abs(dy)/maxaccel) , math.sqrt(6*abs(dz)/Zmaxaccel), 3*abs(dz)/2/0.4)
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
def change_places (dx,dy,dz):
    shutter=0
    shutter_sp=0
    x_p=dx
    x_v=0
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    time=.001+max(math.sqrt(8*abs(dx)/maxaccel) , math.sqrt(8*abs(dy)/maxaccel) , math.sqrt(8*abs(dz)/Zmaxaccel))
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def initialize ():
    return format_pvt([1.23456789,0,0,0,0,0,0,0,0], set_stage, axo)
def write_PVT(PVT):
    #write a PVT file.txt
    f=open('da_SqPill_1w_W'+str(round(pwidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt','w')
    for i in range (0,len(PVT)):
        f.write(" ".join(map(str, PVT[i])))
        f.write('\n')
    f.close()
def make_trajectory():

    PVT=[]
    PVT.append(initialize())
    PVT.append(stage_start_move())
    for q in range (depthcount):
        
        for l in range (solidcount):
            
            for j in range (Ns): #Part 1, solid bar section
                PVT.append(shutter_open(scanv))
                PVT.append(lase (xlength-LW-pwidth-fatchange,scanv))
                PVT.append(shutter_close(scanv))
                PVT.append(flyback1(  (-flydist-2*shutterspace)  /2, -Rshifty/2, Rshiftz/2))
                PVT.append(flyback2(  (-flydist-2*shutterspace)  /2, -Rshifty/2, Rshiftz/2))
                
            for j in range (Np): #Part 2, blinking pillars section
                for i in range (stripecount):
                    PVT.append(shutter_open(scanv))
                    PVT.append(lase (spacing-LW -fatchange,scanv))
                    PVT.append(shutter_close(scanv))
                    PVT.append(lase (pwidth-2*shutterspace+LW +fatchange,scanv))
                PVT.append(flyback1(  (-(spacing+pwidth)*stripecount)  /2, -Rshifty/2, Rshiftz/2))
                PVT.append(flyback2(  (-(spacing+pwidth)*stripecount)  /2, -Rshifty/2, Rshiftz/2))

        for j in range (Ns): #End cap
            PVT.append(shutter_open(scanv))
            PVT.append(lase (xlength-LW-pwidth-fatchange,scanv))
            PVT.append(shutter_close(scanv))
            PVT.append(flyback1(  (-flydist-2*shutterspace) /2, -Rshifty/2, Rshiftz/2))
            PVT.append(flyback2(  (-flydist-2*shutterspace) /2, -Rshifty/2, Rshiftz/2))
            
        PVT.append(left_dropdown((solidcount*N+Ns)*Rshifty,-(solidcount*N+Ns)*Rshiftz+zdrill)) 
    PVT.append(change_places(0,0,-zdrill*depthcount)) #return to 0,0,0
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
    p = pvt_tools.PVT('da_SqPill_1w_W'+str(round(pwidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt')
if set_stage == '5D':
    p = pvt_tools_5D.PVT('da_SqPill_1w_W'+str(round(pwidth*1000))+'_S'+str(round(spacing*1000))+'_os'+str(depthcount)+'_fc'+str(round(fatchange*1000))+'tc'+str(round(tallchange*1000))+'.txt')


p.read()
p.plot()


#plt.plot(p.A[:,1])
plt.plot(p.A[:,0])
plt.plot(p.A[:,1])
plt.plot(p.A[:,2])

plt.figure()
plt.plot(p.V[:,1]) #Y
#plt.plot(p.V[:,3]) #Shutter??
plt.plot(p.V[:,2]) #Z
#Max Velocity x-y = 25, z =0.5
#Max Accel, x-y=300, z=1
