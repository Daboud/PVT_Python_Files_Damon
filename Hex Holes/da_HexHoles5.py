# -*- coding: utf-8 -*-
"""
@author: damon.aboud
"""
import pvt_tools
import pvt_tools_5D
import matplotlib.pyplot as plt
import math

#--------------------------PROCESSING PARAMETERS-----------------------------
xmove=0.0 #Shift stage position after patch to facilitate many experiments
scanv=1.0 #scanning velocity
xlength0=0.6 #width of patch
ylength0=0.3 #height of patch
hexheight=0.100 #width of hexagon (from one face to another! Not corner to corner.)
spacing=0.050 #space between hexagons
depthcount=1 #overscans
Angle=0     /(180/math.pi) #angle of sample. For example 45o would be =math.pi/4. Flat=0
LW=0.0109 #line width
ol0=0.0 #overlap, as a fraction of 1
zdrill=0.030 #Depth of a single scan
set_stage = '5D'
axo=1 #Axis orientation: Set 1 if you want to draw lines along +X and then raster-shift in the -Y direction (like reading a page of a book, right then down). 
#                        Set 2 if you want to draw lines along +Y and then raster-shift in the +X direction (up then right)
#-----------------------------------------------------------------------------

shuttertime=0.001
t60=math.tan(math.pi/3)
s60=math.sin(math.pi/3)
c60=math.cos(math.pi/3)
sa=math.sin(Angle)
ca=math.cos(Angle)
LWa=LW/ca
maxaccel=80.
Zmaxaccel=0.8
zdepth=zdrill*(depthcount-1)

hwidth=hexheight/s60
bwidth=hwidth/2 #width of hex base
UCheight=hexheight+spacing
xgap=spacing/s60 #spacing between hexes, expressed in the X-direction. "Sine of valley"

Rtot=hexheight #total scanning length
Rshift0=LWa*(1-ol0)
N0=Rtot/Rshift0
E=round(N0/2-0.5)
N=E*2+1
Rshift=Rtot/N #position shift for each line
Rshifty=Rshift*ca
Rshiftz=Rshift*sa
OLtrue=1-Rshift/LWa #true overlap percentage
N=int(N) #total shifts per UC

dx=Rshift/t60
UCwidth=2*(hexheight+spacing)*s60

stripecount=round(xlength0/UCwidth) #Number of unit cell repetitions in X
solidcount=round(ylength0/UCheight)*2#Number of unit cell repetitions in Y
solidcount=int(solidcount)
stripecount=int(stripecount)

shutterspace=shuttertime*scanv
zflytime=math.sqrt(6*Rshiftz/Zmaxaccel)
xflytime=3*(scanv**2+math.sqrt(scanv**2-(2/3*maxaccel)*-UCwidth*stripecount)) / maxaccel

flydisto2=UCwidth*stripecount/2
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
    time=max(xflytime,zflytime)
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
    time=max(t2,zflytime)
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
    time=max(t2,zflytime)
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
    time=0.001+max(math.sqrt(6*abs(xp)/maxaccel) , math.sqrt(6*abs(yp)/maxaccel) , math.sqrt(6*abs(zp)/Zmaxaccel) , 3*abs(zp)/2/0.4)
    PVT=[time,x_p,x_v,y_p,y_v,z_p,z_v,shutter,shutter_sp]
    return format_pvt(PVT, set_stage, axo)
def initialize ():
    return format_pvt([1.23456789,0,0,0,0,0,0,0,0], set_stage, axo)
def write_PVT(PVT):
    #write a PVT file.txt
    f=open('da_HexHoles_W'+str(hexheight*1000)+'_S'+str(spacing*1000)+'_os'+str(depthcount)+'.txt','w')
    for i in range (0,len(PVT)):
        f.write(" ".join(map(str, PVT[i])))
        f.write('\n')
    f.close()
def make_trajectory():

    PVT=[]
    PVT.append(initialize())
    #SEE PPT FILE IN GITHUB FOLDER TO UNDERSTAND THE STRUCTURE OF THIS CODE
    
    for q in range (depthcount): #Number of overscans
        PVT.append(stage_start_move())
        for l in range (solidcount):  #Number of unit cell repetitions in Y
            for j in range (N): #Number of raster scan adjustents for a row of hexes
                for i in range (stripecount): #Number of unit cell repetitions in X
                    if (j<N/2): #for the first half (widening)
                        PVT.append(lase ((xgap)+LW,scanv))
                        PVT.append(shutter_open(scanv))
                        PVT.append(lase ((bwidth+2*(j)*dx -2*shutterspace)-LW,scanv))
                        PVT.append(shutter_close(scanv))
                        PVT.append(lase (UCwidth-xgap-(bwidth+2*(j)*dx)-2*shutterspace,scanv))
                    else: #for the second half (reducing width)
                        PVT.append(lase ((xgap)+LW,scanv))
                        PVT.append(shutter_open(scanv))
                        PVT.append(lase ((bwidth+2*(N-j-1)*dx -2*shutterspace)-LW,scanv))
                        PVT.append(shutter_close(scanv))
                        PVT.append(lase (UCwidth-xgap-(bwidth+2*(N-j-1)*dx) -2*shutterspace,scanv))
                if (j<N/2-1): #Widening 
                    PVT.append(flyback1(  (-(UCwidth-2*shutterspace)*stripecount-dx)  /2, -Rshifty/2, Rshiftz/2))
                    PVT.append(flyback2(  (-(UCwidth-2*shutterspace)*stripecount-dx)  /2, -Rshifty/2, Rshiftz/2))
                else: #Reducing
                    PVT.append(flyback1(  (-(UCwidth-2*shutterspace)*stripecount+dx)  /2, -Rshifty/2, Rshiftz/2))
                    PVT.append(flyback2(  (-(UCwidth-2*shutterspace)*stripecount+dx)  /2, -Rshifty/2, Rshiftz/2))
                    
            PVT.append(change_places(-dx,N*Rshifty,-N*Rshiftz)) #Return to start location of that stripe
            mult=((l/2)-math.ceil(l/2))*4+1 #decide to move the next row to the left or right of previous
            PVT.append(change_places(mult*(hexheight+spacing)*s60,-(hexheight+spacing)*c60*ca,(hexheight+spacing)*c60*sa)) 
        
        if math.floor(solidcount/2)==solidcount/2: #return to (0,0)
            PVT.append(change_places(0,solidcount*(hexheight+spacing)*c60*ca,-solidcount*(hexheight+spacing)*c60*sa+zdrill)) 
        else: #return to (0,0)
            PVT.append(change_places(   -(hexheight+spacing)*s60    ,   solidcount*(hexheight+spacing)*c60*ca    ,   -solidcount*(hexheight+spacing)*c60*sa+zdrill   )   )                  

    PVT.append(change_places(0,0,-depthcount*zdrill))  #return Z to (0)

    PVT.append(change_places (xmove,0,0)) #execute xmove

    write_PVT(PVT)

    xsum=sum(row[1] for row in PVT)
    ysum=sum(row[3] for row in PVT)
    zsum=sum(row[5] for row in PVT)
    t=sum(row[0] for row in PVT)/60
    print('Overlap % = ',str(format(OLtrue, '.4f')))
    print('Time (min) = ',str(format(t, '.2f')))
    if t>180:
        print('Time (h) = ',str(format(t/60, '.2f')))
    print('Sum X = ',str(format(xsum, '.5f')))
    print('Sum Y = ',str(format(ysum, '.5f')))
    print('Sum Z = ',str(format(zsum, '.5f')))
   
make_trajectory()
# Check results 

if set_stage == '3D':
    p = pvt_tools.PVT('da_HexHoles_W'+str(hexheight*1000)+'_S'+str(spacing*1000)+'_os'+str(depthcount)+'.txt')
if set_stage == '5D':
    p = pvt_tools_5D.PVT('da_HexHoles_W'+str(hexheight*1000)+'_S'+str(spacing*1000)+'_os'+str(depthcount)+'.txt')

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
