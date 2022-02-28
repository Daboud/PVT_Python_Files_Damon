# -*- coding: utf-8 -*-
"""
@author: damon.aboud
"""

import pvt_tools
import math
import pvt_tools_5D

set_stage = '5D'
outline=0
xmove=0

d=0.300/1.635 #This value represents the slope of the side walls. For example, 0.3/1.635 means that
# we measured a 300um difference between the top and bottom width for a 1.635mm sample. It's essentially
# a tangent function.

 #-----------------------------------------------------
T=0.850 #sample thickness
Dbot=0.100 #intended inlet diameter of cylinder
Cheight=0.340 #height of cylinder
Dtop=Dbot+Cheight*d
Rheight=T-Cheight #thickness of rectangle portion
delta=Rheight*d
xlength0=1.500 #length of long side of the rectangular top portion

xlength=xlength0+delta #Width of rectangle at the cylinder/rectange interface
ylength=Dtop+delta

zdrill=0.125 #depth of one overscan layer
depthcount=round(Rheight/zdrill) #overscans for rectangular hammer head
ol=0.5 #overlap between raster scan lines
scanv=2.0 #scanning translational velocity
zdrill_real=0.001 #This is the value the Z stage actually moves between overscans. I found that we got better
# channels when 'zdrill' above was used to calculate how many overscans are required to drill thru the sample,
# but we actually did not move the position of the Z stage, since the 200mm lens has a super long Rayleigh length

depthcountS=round(Cheight/zdrill)+2 #overscans for cylinder
#next time !!!!!!!!!!! replace depthcountS with the +3, and make depthcount +0 again !!
#--------------------------------------------------------

R=Dtop/2
w0=0.038 #um

Stime=0.001

maxaccel=50.0
Zmaxaccel=0.8
zdepth=zdrill*(depthcount-1)

Rshift=w0*(1-ol)
N=ylength/Rshift/2
N=round(N)
N=int(N)

a=Rshift/2/math.pi##########################
Nt=4
laps=R/Rshift
laps=int(laps)

shutterspace=Stime*scanv


def format_pvt(pvt, stage='3D'):
    if stage == '3D':
        return pvt
    elif stage == '5D':
        return pvt[0:7] + [0.0, 0.0] + pvt[7:9] + [0.0, 0.0]
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
    return format_pvt(PVT, set_stage)
def stage_start_moveS(xv,yv):
    time=4*scanv/maxaccel #time calculated manually, a=8mm/s^2
    shutter=0
    shutter_sp=0
    x_p=-shutterspace
    x_v=xv
    y_p=0
    y_v=yv
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def shutter_openS(xv,yv,zv):
    time=Stime #sec
    shutter=5
    shutter_sp=0
    x_v=xv
    x_p=time*xv
    y_v=yv
    y_p=time*yv
    z_v=zv
    z_p=time*zv
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def shutter_closeS(xv,yv,zv):
    time=Stime #sec
    shutter=-5
    shutter_sp=0
    x_v=xv
    x_p=time*xv
    y_v=yv
    y_p=time*yv
    z_v=zv
    z_p=time*zv
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def shutter_open(v):
    time=Stime #sec
    shutter=5
    shutter_sp=0
    x_v=v
    x_p=time*v
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def shutter_close(v):
    time=Stime #sec
    shutter=-5
    shutter_sp=0
    x_v=v
    x_p=time*v
    y_p=0
    y_v=0
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
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
    return format_pvt(PVT, set_stage)
def right_elbow (dy,dz):
    time=max(2.*scanv/maxaccel , math.sqrt(8*abs(dy)/maxaccel) , math.sqrt(8*abs(dz)/Zmaxaccel))
    math.sqrt(8*abs(dz)/Zmaxaccel)
    shutter=0
    shutter_sp=0   #velocity tend to zero for both the z and y stages:
    x_p=0
    x_v=-scanv
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def left_elbow (dy,dz):
    time=max(2.*scanv/maxaccel , math.sqrt(8*abs(dy)/maxaccel) , math.sqrt(8*abs(dz)/Zmaxaccel))
    shutter=0
    shutter_sp=0   #velocity tend to zero for both the z and y stages:
    x_p=0
    x_v=scanv
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def left_dropdown (dy,dz):
    time=max(4.*scanv/maxaccel , math.sqrt(8*abs(dy)/maxaccel) , math.sqrt(8*abs(dz)/Zmaxaccel))
    shutter=0
    shutter_sp=0   #velocity tend to zero for both the z and y stages:
    x_p=0
    x_v=scanv
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def diam (Nc):
    theta=2*math.pi*Nc/Nt
    Rc=    a    *theta
    shutter=0
    shutter_sp=0
    xpc=Rc*math.cos(2*math.pi/Nt*Nc)
    ypc=Rc*math.sin(2*math.pi/Nt*Nc)
    xpp=(Rc-Rshift/4)*math.cos(2*math.pi/Nt*(Nc-1))
    ypp=(Rc-Rshift/4)*math.sin(2*math.pi/Nt*(Nc-1))    
    if Nc==1:
        xpp=0
        ypp=0
    x_p=xpc-xpp
    x_v=-scanv*math.sin(theta)
    y_p=ypc-ypp
    y_v=scanv*math.cos(theta)
    time=(2*math.pi*Rc/Nt/scanv) #sec
    z_p=0
    z_v=0
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def change_places (dx,dy,dz):
    shutter=0
    shutter_sp=0
    x_p=dx
    x_v=0
    y_p=dy
    y_v=0
    z_p=dz
    z_v=0
    time=max(math.sqrt(8*abs(dx)/maxaccel) , math.sqrt(8*abs(dy)/maxaccel) , math.sqrt(8*abs(dz)/Zmaxaccel))+0.001
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage)
def write_PVT(PVT):
    #write a PVT file.txt
    f=open('Hammer_CyHi'+str(int(Cheight*1000))+'_CyBa'+str(int(Dbot*1000))+'_Rx'+str(int(xlength0*1000))+'_Ry'+str(int(Dtop*1000))+'_Zd'+str(int(zdrill_real*1000))+'_T'+str(T)+'.txt','w')
    for i in range (0,len(PVT)):
        f.write(" ".join(map(str, PVT[i])))
        f.write('\n')
    f.close()
def make_trajectory():

    PVT=[]
    for i in range (depthcount):
        PVT.append(stage_start_move())
        for j in range (N):
            PVT.append(shutter_open(scanv))
            PVT.append(lase (xlength,scanv))
            PVT.append(shutter_close(scanv))
            PVT.append(right_elbow (-Rshift,0))
            PVT.append(shutter_open(-scanv))
            PVT.append(lase (-xlength,-scanv))
            PVT.append(shutter_close(-scanv))
            PVT.append(left_elbow (-Rshift,0))       
        PVT.append(change_places(0,Rshift*2*N,zdrill_real))
        
    PVT.append(change_places(xlength/2,-(2*N-1)*Rshift/2,0))
    for j in range (depthcountS):
        PVT.append(stage_start_moveS(scanv, 0))
        PVT.append(shutter_openS(scanv, 0, 0))
        for i in range (Nt*laps):
        # for i in range (Nt):
            PVT.append(diam(i+1))
        PVT.append(shutter_closeS(0, scanv, 0))
        PVT.append(change_places (-laps*Rshift,-shutterspace,zdrill_real))
    PVT.append(change_places(-xlength/2,ylength/2,-zdrill_real*(depthcount+depthcountS)))
    
    write_PVT(PVT)

    xsum=sum(row[1] for row in PVT)
    ysum=sum(row[3] for row in PVT)
    zsum=sum(row[5] for row in PVT)
    # print('Overlap % = ',str(format(ol, '.4f')))
    print('Time (min) = ',str(format(sum(row[0] for row in PVT)/60, '.2f')))
    print('Sum X = ',str(format(xsum, '.5f')))
    print('Sum Y = ',str(format(ysum, '.5f')))
    print('Sum Z = ',str(format(zsum, '.5f')))
   
    print('Cylinder Diameter (base) = ',str(format(Dbot, '.3f')))
    print('Cylinder Diameter (top) = ',str(format(Dtop, '.3f')))
    print('Rectangle Length (base) = ',str(format(xlength-d*Rheight, '.3f')))
    print('Rectangle width (base) = ',str(format(ylength-d*Rheight, '.3f')))
    print('Rectangle Length (top) = ',str(format(xlength, '.3f')))
    print('Rectangle width (top) = ',str(format(ylength, '.3f')))
    print('Rectangle Depth = ',str(format(depthcount*zdrill, '.3f')))
    print('Cylinder Depth = ',str(format((depthcountS-2)*zdrill, '.3f')), '+',  str(format(2*zdrill, '.3f'))  )
       
make_trajectory()

if set_stage == '3D':
    p = pvt_tools.PVT('Hammer_CyHi'+str(int(Cheight*1000))+'_CyBa'+str(int(Dbot*1000))+'_Rx'+str(int(xlength0*1000))+'_Ry'+str(int(Dtop*1000))+'_Zd'+str(int(zdrill_real*1000))+'_T'+str(T)+'.txt')
if set_stage == '5D':
    p = pvt_tools_5D.PVT('Hammer_CyHi'+str(int(Cheight*1000))+'_CyBa'+str(int(Dbot*1000))+'_Rx'+str(int(xlength0*1000))+'_Ry'+str(int(Dtop*1000))+'_Zd'+str(int(zdrill_real*1000))+'_T'+str(T)+'.txt')

p = pvt_tools.PVT('DJham_CyHi'+str(int(Cheight*1000))+'_CyBa'+str(int(Dbot*1000))+'_Rx'+str(int(xlength0*1000))+'_Ry'+str(int(Dtop*1000))+'_Zd'+str(int(zdrill_real*1000))+'_T'+str(T)+'.txt')
p.read()
p.plot()

