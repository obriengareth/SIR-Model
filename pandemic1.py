###############################################################################
"""
messing with SIR model for Ireland
Gareth O'Brien
"""
###############################################################################
# load libraries
import matplotlib . pyplot as plt
import numpy as np
#import time
import pandas as pd
import imageio  
 #from mpl_toolkits.mplot3d import Axes3D  
#from matplotlib import cm

###############################################################################
# define functions
def xy2ij(x,y,dx,dy):
    i=int (np.floor((x-utm_e1)/dx))
    j=int (np.floor((y-utm_n1)/dy))
    return i,j

import time
def tic():
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " s")
    else:
        print("Toc: start time not set")


###############################################################################
# global variables
'''
global utm_e1
global utm_n1
global utm_e2
global utm_n2
'''
###############################################################################
# input parameters
# scale of the model
no=1;
Nn=2;
Nx=100;
Ny=120;
BigPop = 100;

utm_e1=400000
utm_n1=500000
utm_e2=800000
utm_n2=1000000

# pandemic parameters
beta=0.5
R0=4;
gamma = beta/R0;
transport = 5.85;
mob_neighs = 0.85
no_days = 90

###############################################################################
#  varaibles
dx = (utm_e2-utm_e1)/Nx;
dy = (utm_n2-utm_n1)/Ny;

###############################################################################
# load in population data and convert to regular grid
tic()
data = pd.read_csv('population_data_1km2.csv', usecols=["CENT_EAST","CENT_NORTH","TOT_P"] );#,names=['x','y','z'])
Ld, nr =data.shape

X = np.zeros((Nx,Ny), dtype=np.float32)
Y = np.zeros((Nx,Ny), dtype=np.float32)
Pij = np.zeros((Nx*Ny,2), dtype=np.float32)

pop = np.zeros((Nx,Ny), dtype=np.float32)
for k in range(0,Ld):
    i, j = xy2ij(data.CENT_EAST[k],data.CENT_NORTH[k],dx,dy)
    pop[i][j] = pop[i][j] + data.TOT_P[k];
toc()
############################################################################### 
# build coastline for plotting
coast = np.zeros((Ny*2,2), dtype=np.float32)
mask = np.zeros((Nx,Ny), dtype=np.float32)
mask[:][:]=1
k=0;
for j in range(0,Ny):
    for i in range(1,Nx-1):
        if(pop[i][j]!=0 and pop[i-1][j]==0 and pop[i+1][j]!=0):
            coast[k][0]=i*dx+utm_e1
            coast[k][1]=j*dy+utm_n1
            k=k+1;
            break
for j in range(Ny-1,0,-1):
    for i in range(Nx-1,1,-1):
        if(pop[i][j]!=0 and pop[i-1][j]!=0 and pop[i+1][j]==0):
            coast[k][0]=i*dx+utm_e1
            coast[k][1]=j*dy+utm_n1
            k=k+1;
            break
Lco=k;   

'''
fig = plt.figure(num=no,figsize=(9, 6)) ;
plt.clf();
plt.pcolor(X,Y,pop,vmin=0,vmax=1);
plt.plot(coast[0:Lco,0],coast[0:Lco,1],'r');
plt.xlabel('UTM X')
plt.ylabel('UTM Y')
#plt.title(r'Pore Pressure Evolution')
plt.pause(0.1);
plt.show() 
 '''    
############################################################################### 
# build lists of population and define interaction cells by labelling indices
P = np.zeros((Nx*Ny,1), dtype=np.float32)      
k=0;
for i in range(0, Nx):
    for j in range(0,Ny):
        X[i][j] = i*dx + utm_e1;
        Y[i][j] = j*dy + utm_n1;   
        P[k][0] = pop[i][j];
        Pij[k][0]=i;
        Pij[k][1]=j;
        k = k +1 ;
Lp=k;

Pindex = np.argsort(P,axis=0)

AllNeigh =   (2*Nn+1)*(2*Nn+1) + BigPop ; 
Neigh = np.zeros((Nx*Ny,AllNeigh), dtype=np.float32)
for k in range(0,Lp):
    l=0;
    for i in range(-Nn,Nn+1):
        for j in range(-Nn,Nn+1):
            Neigh[k,l]=(Pij[k][0]-i)*Ny + (Pij[k][1]-j);
            if(Neigh[k][l]>=Lp):
                Neigh[k][l]=0;
            l=l+1;
            
    for m in range(0, BigPop):
        Neigh[k][l+m]  = Pindex[Lp-m-1][0] ;

###############################################################################
# build a normalised mobility network
mob = np.zeros((Nx*Ny,AllNeigh), dtype=np.float32)
for k in range(0,Lp):
    l=0;
    for m in range(0,(2*Nn+1)*(2*Nn+1)):
        mm=int (Neigh[k][m]);
        if(k==mm):
            mob[k][m]  = 0 ;        
        else:
            mob[k][m]  = mob_neighs ;        

    for m in range((2*Nn+1)*(2*Nn+1),AllNeigh):
        mm=int (Neigh[k][m]);
        if(k==mm):
            mob[k][m]  = 0 ;        
        else:
            mob[k][m]  = ( P[k][0]*P[mm][0] )/ (1.0*P[Pindex[Nx*Ny-1][0]][0]) / (1.0*P[Pindex[Nx*Ny-1][0]][0]) ;   

  
###############################################################################
# doing time loop over the days
S = np.zeros((Lp,1), dtype=np.float32)
I = np.zeros((Lp,1), dtype=np.float32)
R = np.zeros((Lp,1), dtype=np.float32)
S[:,0]=P[:,0];
np.sum(S)
dummy = np.zeros((Nx,Ny), dtype=np.float32)
tt = np.zeros((no_days,), dtype=np.float32)
Ss = np.zeros((no_days,), dtype=np.float32)
Is = np.zeros((no_days,), dtype=np.float32)
Rs = np.zeros((no_days,), dtype=np.float32)

for time_days in range(0,no_days):
    tt[time_days]=time_days

for m in range(Lp-20,Lp):
    I[Pindex[m][0]][0] = round(P[Pindex[m][0]][0]*0.00007);
    S[Pindex[m][0]][0] = P[Pindex[m][0]][0] - I[Pindex[m][0]][0];

fig = plt.figure(num=no,figsize=(12, 6)) ;
plt.clf(); 
gs = plt.GridSpec(3, 3,hspace=0.5, wspace=0.3);#width_ratios=[1, 3], height_ratios=[2, 1])
ax1 = fig.add_subplot(gs[0:2, 0:1])
ax2 = fig.add_subplot(gs[0:2 ,1:2])
ax3 = fig.add_subplot(gs[2:3 ,0:3])
ax4 = fig.add_subplot(gs[0:2 ,2:3])
images=[]


for time_days in range(0,no_days):
    for k in range(0,Lp):
        t1=0; t2=0; t3=0; t4=0;
        
        for m in range(0,AllNeigh):
            mm = int(Neigh[k][m]);
            if(P[mm][0]!=0):
                t1 = t1 + transport*S[k][0]*mob[k][m]*I[mm][0]*beta/P[mm][0]
                t2 = t2 + mob[k][m]
        if( (P[k][0]+t2)==0 ):
            t3=0;
        else:
            t3 = t1/(P[k][0]+t2);
    
        I0=I[k][0]; 
        if(P[k][0]==0):
            t4=0;
        else:
            t4=beta*S[k][0]*I0/P[k][0] ;
  
        R[k][0] = R[k][0] + gamma*I0;
        I[k][0] = I[k][0] + t4 + t3 - gamma*I0;
        S[k][0] = S[k][0] - t4 - t3;

    print('Day ',time_days,'I= ',np.sum(I),' R= ',np.sum(R))
    ################################################################
    # plotting figure
    fig = plt.figure(num=no,figsize=(12, 6)) ;
    plt.clf(); 
    gs = plt.GridSpec(3, 3,hspace=0.5, wspace=0.3);#width_ratios=[1, 3], height_ratios=[2, 1])
    ax1 = fig.add_subplot(gs[0:2, 0:1])
    ax2 = fig.add_subplot(gs[0:2 ,1:2])
    ax3 = fig.add_subplot(gs[2:3 ,0:3])
    ax4 = fig.add_subplot(gs[0:2 ,2:3])


    for k in range(0,Lp):  
        i=int(Pij[k][0]);
        j=int(Pij[k][1]);
        dummy[i][j] = P[k][0];
                
    pcm=ax1.pcolor(X,Y,dummy,vmin=0,vmax=1e3,cmap='Purples');
    fig.colorbar(pcm,ax=ax1);
    ax1.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
    ax1.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
    ax1.set_ylabel('North')
    ax1.set_xlabel('East')
    ax1.grid()
    ax1.set_title("Population")


    for k in range(0,Lp):  
        i=int(Pij[k][0]);
        j=int(Pij[k][1]);
        dummy[i][j] = I[k][0];
                
    pcm=ax2.pcolor(X,Y,dummy,vmin=0,vmax=1e3,cmap='Reds');
    fig.colorbar(pcm,ax=ax2);
    ax2.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
    #ax2.set_ylabel('North')
    #ax2.set_xlabel('East')
    #ax2.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.grid()
    buf = "# Infected on Day %d " % time_days
    ax2.set_title(buf)
    

    for k in range(0,Lp):  
        i=int(Pij[k][0]);
        j=int(Pij[k][1]);
        dummy[i][j] = R[k][0];
                
    pcm=ax4.pcolor(X,Y,dummy,vmin=0,vmax=1e3,cmap='Greens');
    fig.colorbar(pcm,ax=ax4);
    ax4.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
    #ax4.set_ylabel('North')
    ax4.set_xlabel('East')
    ax4.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
    ax4.set_yticks([])
    #ax4.grid()
    buf = "# Recovered on Day %d " % time_days
    ax4.set_title(buf)
    
    Ss[time_days] = np.sum(S)
    Is[time_days] = np.sum(I)
    Rs[time_days] = np.sum(R)
    
    ax3.plot(tt[0:time_days],Is[0:time_days],'r');    
    ax3.plot(tt[0:time_days],Rs[0:time_days],'g');
    ax3.plot(tt[0:time_days],Ss[0:time_days],'b');
    ax3.set_ylabel('Population')
    ax3.set_xlabel('Time (Days)')
    ax3.set_title("SIR Curves")
    ax3.axis([0,90,0,5e6])
    ax3.grid()
    
    plt.pause(0.2);
    plt.show()
    plt.savefig('mess.png',dpi=120)            
    images.append(imageio.imread('mess.png'))


imageio.mimsave('Ireland_CORVID_M2.gif', images,'GIF',duration=0.2)
    
###############################################################################
# plot test figures        

'''
fig = plt.figure(num=no,figsize=(12, 6)) ;
plt.clf(); 
gs = plt.GridSpec(3, 3,hspace=0.5, wspace=0.3);#width_ratios=[1, 3], height_ratios=[2, 1])
ax1 = fig.add_subplot(gs[0:2, 0:1])
ax2 = fig.add_subplot(gs[0:2 ,1:2])
ax3 = fig.add_subplot(gs[2:3 ,0:3])
ax4 = fig.add_subplot(gs[0:2 ,2:3])

for k in range(0,Lp):  
    i=int(Pij[k][0]);
    j=int(Pij[k][1]);
    dummy[i][j] = P[k][0];
            
pcm=ax1.pcolor(X,Y,dummy,vmin=0,vmax=1e3,cmap='Purples');
fig.colorbar(pcm,ax=ax1);
ax1.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
ax1.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
ax1.set_ylabel('North')
ax1.set_xlabel('East')
ax1.grid()
ax1.set_title("Population")


for k in range(0,Lp):  
    i=int(Pij[k][0]);
    j=int(Pij[k][1]);
    dummy[i][j] = I[k][0];
            
pcm=ax2.pcolor(X,Y,dummy,vmin=0,vmax=1e3,cmap='Reds');
fig.colorbar(pcm,ax=ax2);
ax2.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
#ax2.set_ylabel('North')
ax2.set_xlabel('East')
#ax2.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
ax2.set_yticks([])
ax2.grid()
buf = "# Infected on Day %d " % time_days
ax2.set_title(buf)


for k in range(0,Lp):  
    i=int(Pij[k][0]);
    j=int(Pij[k][1]);
    dummy[i][j] = R[k][0];
            
pcm=ax4.pcolor(X,Y,dummy,vmin=0,vmax=1e4,cmap='Greens');
fig.colorbar(pcm,ax=ax4);
ax4.plot(coast[0:Lco,0],coast[0:Lco,1],'k');
#ax4.set_ylabel('North')
ax4.set_xlabel('East')
ax4.set_xticks([utm_e1,(utm_e1+utm_e2)*0.5, utm_e2])
ax4.set_yticks([])
ax4.grid()
buf = "# Recovered on Day %d " % time_days
ax4.set_title(buf)

Ss[time_days] = np.sum(S)
Is[time_days] = np.sum(I)
Rs[time_days] = np.sum(R)

ax3.plot(tt[0:time_days],Is[0:time_days],'r');    
ax3.plot(tt[0:time_days],Rs[0:time_days],'g');
ax3.plot(tt[0:time_days],Ss[0:time_days],'b');
ax3.set_ylabel('Population')
ax3.set_xlabel('Time (Days)')
ax3.set_title("SIR Curves")
ax3.axis([0,90,0,1e6])
ax3.grid()
plt.pause(0.1);
plt.show()

'''