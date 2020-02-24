import matplotlib.pyplot as plt


import numpy as numpy
from scipy.integrate import odeint
plt.interactive(True)
from matplotlib.widgets import Slider, Button
from math import*

class Bob():

   
    def __init__(self,m,inang,l):
        #in puts are mass, length and initial angle
        #length
        self.l = l
        #mass
        self.m = m
        #angular displacement
        self.theta = inang
        #angular velocity
        self.v = 0
        #displacement in the x direction
        self.x = 0
        #displacement in the y direction
        self.y = 0
        #potential energy
        self.p = 0
        #mechanical energy
        self.energy = 0
        #kinetic energy
        self.ke = 0
        #potential energy
        self.pe = 0
       
       




def mainloop(l1,l2,m1,m2,inang1,inang2,tend,dt):
    """
       Runs the 'main loop'
   
    """
   
    p1 = Bob(l1,m1,inang1)
    p2 = Bob(l2,m2,inang2)
    p1.a,p2.a = getacc(p1,p2)
   
   
    fig,axes,data = initialize_plots(p1,p2)
   
    t = numpy.arange(int(tend/dt)+1)*dt


    i = 1
    p1t=slidert1.val
    p2t=slidert2.val
    ylimholder02=0
    ylimholder12=0
    c102=0
   
   
   
    for i,ti in enumerate(t[1:],start=1):
        c102+=1

        p1.a,p2.a = getacc(p1,p2)
        lsv=[p1.v,p2.v]
        lsv.sort()
        ls=[p1.a,p2.a,-p1.a,-p2.a]
        ls.sort()

        p1.ma=sliderm1.val
        p2.m=sliderm2.val
        p1.l=sliderl1.val
        p2.l=sliderl2.val
       
       
        if (p1t==slidert1.val)==False:
            p1.theta=slidert1.val
            p1t=slidert1.val
        if (p2t==slidert2.val)==False:
            p2.theta=slidert2.val
            p2t=slidert2.val
       
        #g=sliderg.val
        axes[0,0].set_ylim(-(p2.l+p1.l),(p1.l+p2.l))
        axes[0,0].set_xlim(-(p2.l+p1.l),(p1.l+p2.l))
        axes[1,1].set_ylim(-(p1.energy+p2.energy)*1.25,(p1.energy+p2.energy)*1.25)
        if c102%6==0:
            axes[1,2].set_ylim(ls[0]*1.25,ls[3]*1.25)
        if c102%1==0:
            if p1.a<0 and p2.a<0:
                exit
            axes[0,2].set_ylim(-abs((lsv[1]+lsv[0])*1.25),abs(lsv[1]+lsv[0])*1.25)
       
        #while p1.a>ylimholder02 or p2.a>ylimholder02 or p1.a<-ylimholder02 or p2.a<-ylimholder02:
        #  
        #    axes[1,2].set_ylim(-(c102),(c102))
        #    
        #    c102+=0.1
           
        #if p1.a>p2.a:
        #    axes[1,2].set_ylim(-(p1.energy+p2.energy)*1.25,(p1.energy+p2.energy)*1.25)
       
        integrator(p1,p2,dt)
        x1,y1,x2,y2 = getpos(p1,p2)
        e1,e2 = get_energies(p1,p2)
        update_plots(p1,p2,ti,fig,axes,data)
   

    #return testslider




def getacc(p1,p2):
    ''' Just the formula to get this , to see derivation , please refer to notes'''
    m1=p1.m
    m2=p2.m
    r1=p1.l
    r2=p2.l
    theta1_dot=p1.v
    theta2_dot=p2.v
    theta1=p1.theta
    theta2=p2.theta

    num1 = (-g*(2*m1+m2)*sin(theta1))
    num2 = -m2*g*sin(theta1-2*theta2)
    num3 = (-2*sin(theta1-theta2)*m2*(theta2_dot**2*r2+theta1_dot**2*r1*cos(theta1-theta2)))
    denom1 = r1*(2*m1+m2-m2*cos(2*theta1-2*theta2))
    theta1_dotdot = (num1 + num2 + num3)/denom1

    num4 = 2*sin(theta1-theta2)
    num5 = (theta1_dot**2*r1*(m1+m2))
    num6 = g*(m1+m2)*cos(theta1)
    num7 = theta2_dot**2*r2*m2*cos(theta1-theta2)
    denom2 = r2*(2*m1+m2-m2*cos(2*theta1-2*theta2))
    theta2_dotdot = (num4*(num5+num6+num7))/denom2

    return theta1_dotdot , theta2_dotdot


   

def getpos(p1,p2):
    """ Calculate the cartesian displacement per bob """

    x1,x2 = getposx(p1,p2)
    y1,y2 = getposy(p1,p2)

    p1.x = x1
    p1.y = y1
    p2.x = x2
    p2.y = y2
   

   

    return x1,y1,x2,y2

def getposx(p1,p2):
    """ Calculate the cartesian displacement per bob in x"""
    l2 = p2.l
    t2 = p2.theta
   
    l1 = p1.l
    t1 = p1.theta
   
    x1 = l1*numpy.sin(t1)
    x2 = x1 + l2*numpy.sin(t2)
   
    return x1 , x2
   
def getposy(p1,p2):
    """ Calculate the cartesian displacement per bob in y"""
 
    l1 = p1.l
    t1 = p1.theta

    l2 = p2.l
    t2 = p2.theta
   
    y1 = -l1*numpy.cos(t1)
    y2 = y1 - l2*numpy.cos(t2)
    return y1 , y2

   


def get_energies(p1,p2):
    """***useful even if not in lagranian***energies for axes[1,1]"""
    x1,y1,x2,y2 = getpos(p1,p2)

    vx1 = -y1*p1.v
    vy1 = x1*p1.v
    vx2 = vx1 + (y1-y2)*p2.v
    vy2 = vy1 + (x2-x1)*p2.v

    p1.ke = 0.5*p1.m*(vx1**2 + vy1**2)
    p1.pe = p1.m*g*y1
    p1.energy = p1.ke + p1.pe
    p2.ke = .5*p2.m*(vx2**2 + vy2**2)
    p2.pe = p2.m*g*y2
    p2.energy = p2.ke + p2.pe
    return p1.energy,p2.energy

def shortgetacc(integratelist,t,p1,p2):
    """ gets newtonian acceleration, but compiles it , also"""

    p1.theta,p2.theta,p1.v,p2.v = integratelist
    a1,a2 = getacc(p1,p2)
    paslis = numpy.array([p1.v,p2.v,a1,a2])
    return paslis


def integrator(p1,p2,dt):
    """ Integation and bounds"""
   
    integratelist = numpy.array([p1.theta,p2.theta,p1.v,p2.v])

    paslis=odeint(shortgetacc,integratelist,[0,dt],args=(p1,p2))
    p1.theta,p2.theta,p1.v,p2.v = paslis[1]

           
    kpsngwthnbnds1(p1)
    kpsngwthnbnds2(p2)
    return


def kpsngwthnbnds1(p1):
    if p1.theta > numpy.pi:
        while p1.theta > numpy.pi:
            p1.theta -= 2*(numpy.pi)
           
    if p1.theta < -numpy.pi:
        while p1.theta < -numpy.pi:
            p1.theta += (2*numpy.pi)
    return

def kpsngwthnbnds2(p2):
    if p2.theta > numpy.pi:
        while p2.theta > numpy.pi:
            p2.theta -= (2*numpy.pi)
           
    if p2.theta < -numpy.pi:
        while p2.theta < -numpy.pi:
            p2.theta += (2*numpy.pi)
    return

def initialize_plots(p1,p2):
    """ Set the plots that will be used by the animation """
       
   
    fig,axes = plt.subplots(2,3,figsize=(12,8))
    # Figure with four subplots - matplotlib - don't add subplots/axes to fig
   
   
   

    xlist = [0,p1.x,p2.x]      
    ylist = [0,p1.y,p2.y]
    data =[]

    axes[0,0].plot([-0.5,0.5],[0,0],'-k',linewidth=2)
   
   

    line, = axes[0,0].plot(xlist,ylist, '-bo', markersize=7,linewidth=1)
   
    line1, = axes[0,0].plot(p1.x,p1.y,'-g',linewidth=1)
    line2, = axes[0,0].plot(p2.x,p2.y,'-r',linewidth=1)
    axes[0,0].set_xlim(-(p1.l + p2.l),p1.l + p2.l)
   
   

    axes[0,0].set_ylim(-(p1.l + p2.l),p1.l + p2.l)
    axes[0,0].set_title('t = 0',fontsize=15)
    axes[0,0].set_xlabel('x (m)',fontsize=15)
    axes[0,0].set_ylabel('y (m)',fontsize=15)
    data.append([line,line1,line2])
   
    line2, = axes[1,0].plot(0,p1.theta,'-g',label='bob1')
    line3, = axes[1,0].plot(0,p2.theta,'-r',label='bob2')

    axes[1,0].set_ylim(-numpy.pi,numpy.pi)
    axes[1,0].set_xlabel('t(s)',fontsize=15)
    axes[1,0].set_ylabel('$\\theta(rad)$',fontsize=15)
    data.append([line2,line3])
    axes[1,0].legend()

   
   



    line1, = axes[0,1].plot(p1.theta,p1.v,'b.',label='bob1')

    axes[0,1].set_xlabel('$\\theta(rad)$',fontsize=15)
    axes[0,1].set_ylabel('$\\omega(rad/s)$',fontsize=15)
    axes[0,1].set_xlim(-4.14,4.14)
    axes[0,1].set_ylim(-6,6)
   


    line2, = axes[0,1].plot(p2.theta,p2.v,'r.',label='bob2')
    axes[0,1].legend()

    axes[1,1].set_xlabel('t (s)',fontsize=15)
    axes[1,1].set_ylabel('Energy(J)',fontsize=15)

    data.append([line1,line2])
   
    line1, = axes[1,1].plot(0,p1.energy,'-g',label='bob1')
    line2, = axes[1,1].plot(0,p2.energy,'-r',label='bob2')
    line3, = axes[1,1].plot(0,p1.energy+p2.energy,'-m',label='total')
    data.append([line1,line2,line3])
    axes[1,1].legend()
   

    global slideraxesm1
    global sliderm1
    slideraxesm1=plt.axes([0.04,0.93,0.07,0.05])
    sliderm1=Slider(slideraxesm1,"mass1",valmin=0,valmax=10)

    global slideraxesm2
    global sliderm2
    slideraxesm2=plt.axes([0.18,0.93,0.07,0.05])
    sliderm2=Slider(slideraxesm2,"mass2",valmin=0,valmax=10)

    global slideraxesl1
    global sliderl1
    slideraxesl1=plt.axes([0.36,0.93,0.07,0.05])
    sliderl1=Slider(slideraxesl1,"length1",valmin=0.01,valmax=1)

    global slideraxesl2
    global sliderl2
    slideraxesl2=plt.axes([0.54,0.93,0.07,0.05])
    sliderl2=Slider(slideraxesl2,"length2",valmin=0.01,valmax=1)

   
    global slideraxest1
    global slidert1
    slideraxest1=plt.axes([0.72,0.93,0.07,0.05])
    slidert1=Slider(slideraxest1,"theta1",valmin=0,valmax=3.14)

    global slideraxest2
    global slidert2
    slideraxest2=plt.axes([0.90,0.93,0.07,0.05])
    slidert2=Slider(slideraxest2,"theta2",valmin=0,valmax=3.14)

   
    axes[1,2].set_xlabel('t(s)',fontsize=12)
    axes[1,2].set_ylabel('a (ms^-2)',fontsize=6)
    line2, = axes[1,2].plot(0,p1.theta,'-g',label='bob1')
    line3, = axes[1,2].plot(0,p2.theta,'-r',label='bob2')
    axes[1,2].set_ylim(-numpy.pi,numpy.pi)
    axes[1,2].legend()
    data.append([line2,line3])#data insert 4
   
   
    axes[0,2].set_xlabel('t(s)',fontsize=12)
    axes[0,2].set_ylabel('$\\omega(rad/s)$',fontsize=6)
    line2, = axes[0,2].plot(0,p1.theta,'-g',label='bob1')
    line3, = axes[0,2].plot(0,p2.theta,'-r',label='bob2')
    axes[0,2].legend()
    axes[0,2].set_ylim(-numpy.pi,numpy.pi)
    data.append([line2,line3])
   

   

   
    plt.show()
   
   
    return fig,axes,data







def update_plots(p1,p2,t,fig,axes,data):
    """ Updates all of the plots. """
    tstr=str(t)
    c=0
    holder=("The time is -->")
    for i in tstr:
        c+=1
        if c<=4:
            holder+=i
   
   
    axes[0,0].set_title(holder)
    line,line1,line2 = data[0]
    line.set_xdata([0,p1.x,p2.x])
    line.set_ydata([0,p1.y,p2.y])
    line1.set_xdata( numpy.append(line1.get_xdata(),p1.x))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1.y))
    line2.set_xdata( numpy.append(line2.get_xdata(),p2.x))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2.y))

    line1,line2 = data[1]
    line1.set_xdata( numpy.append(line1.get_xdata(), t))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1.theta))
    line2.set_xdata( numpy.append(line2.get_xdata(), t))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2.theta))
    if t > axes[1,0].get_xlim()[1]:
        axes[1,0].set_xlim(0,2*t)



    line1,line2 = data[2]
    line1.set_xdata( numpy.append(line1.get_xdata(), p1.theta))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1.v))
    line2.set_xdata( numpy.append(line2.get_xdata(), p2.theta))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2.v))


    line1,line2,line3 = data[3]
    line1.set_xdata( numpy.append(line1.get_xdata(), t))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1.energy))
    line2.set_xdata( numpy.append(line2.get_xdata(), t))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2.energy))
    line3.set_xdata( numpy.append(line3.get_xdata(), t))
    line3.set_ydata(numpy.append(line3.get_ydata(),p1.energy+p2.energy))
    if t > axes[1,1].get_xlim()[1]:
        axes[1,1].set_xlim(0,2*t)

    ######test######
    p1acc , p2acc = getacc(p1 , p2)
    line1,line2 = data[4]
    line1.set_xdata( numpy.append(line1.get_xdata(), t))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1acc))
    line2.set_xdata( numpy.append(line2.get_xdata(), t))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2acc))
    if t > axes[0,2].get_xlim()[1]:
        axes[0,2].set_xlim(0,2*t)
    ###
    line1,line2 = data[5]
    line1.set_xdata( numpy.append(line1.get_xdata(), t))
    line1.set_ydata(numpy.append(line1.get_ydata(),p1.v))
    line2.set_xdata( numpy.append(line2.get_xdata(), t))
    line2.set_ydata(numpy.append(line2.get_ydata(),p2.v))
    if t > axes[1,2].get_xlim()[1]:
        axes[1,2].set_xlim(0,2*t)


    plt.pause(1*10**(-7))
    plt.show()
    return

global g
g=9.81
a=mainloop( 0.5 , 0.5 , 0.5 , 0.5 , 1, 1, 20, 0.005)