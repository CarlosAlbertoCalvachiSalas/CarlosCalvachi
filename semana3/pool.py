import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as anim
from tqdm import tqdm 
from scipy.optimize import curve_fit


class Particle():
    
    # init
    def __init__(self, r0,v0,a0,t,m,radius,Id):
        
        self.dt  = t[1] - t[0]
        
        self.r = r0
        self.v = v0
        self.a = a0
        
        self.rVector = np.zeros( (len(t),len(r0)) )
        self.vVector = np.zeros( (len(t),len(v0)) )
        self.aVector = np.zeros( (len(t),len(a0)) )
        
        self.m = m
        self.radius = radius
        self.Id = Id
        
    # Method
    def Evolution(self,i):
        
        self.SetPosition(i,self.r)
        self.SetVelocity(i,self.v)
        self.SetAcceleration(i, self.a)
        
        # Euler method
        self.r += self.dt * self.v
        self.v += self.dt * self.a
    
    def CheckWallLimits(self,limits,dim=2, restitution = 1, method = 0):

        if(method == 1): # CheckWallLimits de la clase
            for i in range(dim):

                if (self.r[i] + self.radius  > limits[i]): 
                    self.v[i] = - restitution*self.v[i]

                if (self.r[i] - self.radius < - limits[i]):
                    self.v[i] = - restitution*self.v[i]
         
        else: # CheckWallLimits alternativo
            for i in range(dim):

                if (self.r[i] + self.radius  > limits[i] and self.v[i] > 0): 
                    self.v[i] = - restitution*self.v[i]

                if (self.r[i] - self.radius < - limits[i] and self.v[i] < 0):
                    self.v[i] = - restitution*self.v[i]
            
        
        return (self.r[i] + self.radius  > limits[i] and self.v[i] > 0) or (self.r[i] - self.radius < - limits[i] and self.v[i] < 0)
    # Setters
    def SetPosition(self,i,r):
        self.rVector[i] = r
        
    def SetVelocity(self,i,v):
        self.vVector[i] = v

    def SetAcceleration(self,i,a):
        self.aVector[i] = a

    # Getters  
    def GetPositionVector(self):
        return self.rVector
    
    def GetRPositionVector(self):
        return self.RrVector 

    def GetVelocityVector(self):
        return self.vVector

    def GetRVelocityVector(self):
        return self.RvVector 
    
    def GetKineticEnergy(self):
        return (1/2)*self.m*np.sum(self.vVector ** 2, axis = 1 )

    def GetPotentialEnergy(self):
        return self.m*9.8*self.rVector[:, 1]

    def GetMechanicalEnergy(self):
        return self.GetPotentialEnergy() + self.GetKineticEnergy()

    def GetR(self):
        return self.radius
    
    def ReduceSize(self,factor):
        
        self.RrVector = np.array([self.rVector[0]]) # initial condition
        self.RvVector = np.array([self.vVector[0]]) # initial condition
        
        for i in range(1,len(self.rVector)):
            if i%factor == 0:
                self.RrVector = np.vstack([self.RrVector,self.rVector[i]])
                self.RvVector = np.vstack([self.RvVector,self.vVector[i]])
                
               # print(self.RrVector)

def GetParticles(NParticles,Limit,Velo,Dim=2,dt=0.1):
    
    Particles_ = []
    
    for i in range(NParticles):
        
        x0 = np.random.uniform( -Limit+1.0, Limit-1.0, size=Dim )
        v0 = np.random.uniform( -Velo, Velo, size=Dim)
        a0 = np.zeros(Dim)
        
        p = Particle(x0,v0,a0,t,1.,1.0,i)
        
        Particles_.append(p)
        
    return Particles_


def contactForceParticles(p1, p2, K = 100):
    x1 = np.array(p1.r)
    x2 = np.array(p2.r)

    rvector = x1 - x2

    rnorm = np.sqrt(np.sum(rvector**2))

    r3vector = rvector * (rnorm**2)

    Fmagnitude = K*r3vector*(rnorm < p1.GetR()+p2.GetR()) 

    F21 = -Fmagnitude
    F12 = Fmagnitude

    p1.a = F12/p1.m
    p2.a = F21/p2.m

    return p1.a

def energyPerParticle(p1, p2, K = 100):

    x1 = np.array(p1.r)
    x2 = np.array(p2.r)

    rvector = x1 - x2

    rnorm = np.sqrt(np.sum(rvector**2))

    return (1/8)*K*(rnorm**4)*(rnorm < p1.GetR()+p2.GetR()) 

def RunInteractiveSimulation(t, Particles, Dim = 2):

    potentialEnergy = np.zeros(len(t))
    kineticEnergy   = np.zeros(len(t))

    for it in tqdm(range(len(t))): # Evolucion temporal
        
        accelerationArray = np.zeros([len(Particles), Dim])

        for i in range(len(Particles)):
            otherIndexesArray = np.arange(len(Particles))
            otherIndexesArray = np.append(otherIndexesArray[:i], otherIndexesArray[i+1:])

            for j in otherIndexesArray:
                accelerationArray[i] = accelerationArray[i] + contactForceParticles(Particles[i], Particles[j])
                potentialEnergy[it] = potentialEnergy[it] + energyPerParticle(Particles[i], Particles[j])

        for i in range(len(Particles)):
            Particles[i].a = accelerationArray[i]
            kineticEnergy[it] = kineticEnergy[it] + (1/2)*Particles[i].m*np.sum(Particles[i].v**2)
            Particles[i].CheckWallLimits(Limits, dim = Dim)
            Particles[i].Evolution(it)
        
    totalEnergy = potentialEnergy + kineticEnergy

    plt.plot(t, totalEnergy, label = 'Energía Total', linestyle = '--')
    plt.plot(t, kineticEnergy, label = 'Energía Cinética', linestyle = '--')
    plt.plot(t, potentialEnergy, label = 'Energía Potencial', linestyle = '--')
    plt.legend()
    plt.show()
    
    return Particles

  
def RunSimulation(t,NParticles = 100, Velo = 6, Dim = 2):
    
    Particles = GetParticles(NParticles,Limits[0],Velo = Velo, dt=dt, Dim = Dim)
    
    for it in tqdm(range(len(t))): # Evolucion temporal
        for i in range(len(Particles)):
            Particles[i].CheckWallLimits(Limits, dim = Dim)
            Particles[i].Evolution(it)
       
    return Particles

def RunDampedBallSimulation(t, r0, v0, a0, e, Limits, tmax = 30, method = 1):
    p = Particle(r0, v0, a0, t, 1, 1, 0)

    for tStep in tqdm(range(len(t))):
        p.CheckWallLimits(Limits, restitution=e, method = method)
        p.Evolution(tStep)

    return p

def ReduceTime(t,factor, Particles = []):
    
    for p in Particles:
        p.ReduceSize(factor)
        
    Newt = []
    
    for i in range(len(t)):
        if i%factor == 0:
            Newt.append(t[i])
            
    return np.array(Newt)


Limits = np.array([20.,20.])
dt = 1E-4
tmax = 10
t = np.arange(0,tmax+dt,dt)


p1 = Particle(r0 = np.array([-10., 1.5 ]), v0 =  np.array([20, 0. ]),a0 =  np.array([0. ,0. ]), t = t , m = 1,radius = 2, Id = 0)
p2 = Particle(r0 = np.array([0., -1.6 ]), v0 =  np.array([0., 0. ]),a0 =  np.array([0. ,0. ]), t = t , m = 1,radius = 2, Id = 1)
p3 = Particle(r0 = np.array([-15., -15. ]), v0 =  np.array([0., 0. ]),a0 =  np.array([0. ,0. ]), t = t , m = 1,radius = 2, Id = 2)

Particles = np.array([p1, p2, p3])

simulatedParticles = RunInteractiveSimulation(t, Particles)

redt = ReduceTime(t, 500, Particles = Particles)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1)

def init():
    ax.set_xlim(-Limits[0],Limits[0])
    ax.set_ylim(-Limits[1],Limits[1])
    ax.set_xlabel('x(m)')
    ax.set_ylabel('y(m)')

colors = ['olive', 'red', 'black']
def Update(i, Particles = []):
    
    plot = ax.clear()
    init()
    ax.set_title(r'$t=%.2f \ seconds$ CheckWallLimits ' %(redt[i]), fontsize=15)
    
    for p in Particles:
        x = p.GetRPositionVector()[i,0]
        y = p.GetRPositionVector()[i,1]
     
        vx = p.GetRVelocityVector()[i,0]
        vy = p.GetRVelocityVector()[i,1]
        
        circle = plt.Circle((x,y), p.GetR(), color=colors[p.Id], fill=True, label = p.Id)
        ax.add_patch(circle)
        ax.arrow(x,y,vx,vy,color='r',head_width=0.5)
        plt.legend(loc = 'upper right')
        
    return plot

Animation = anim.FuncAnimation(fig, lambda i: Update(i, Particles = Particles), frames=len(redt), init_func=init, interval = 10, repeat = False)

plt.show()
plt.close()





"""r0 = np.array([-15., 5.])
v0 = np.array([1., 0.])
a0  = np.array([0., -9.8])

dampedParticle = RunDampedBallSimulation(t, r0, v0, a0, 0.9, Limits, 30, method = 1)

r0 = np.array([-15., 5.])
v0 = np.array([1., 0.])
a0  = np.array([0., -9.8])

dampedParticle2 = RunDampedBallSimulation(t, r0, v0, a0, 0.9, Limits, 30, method = 0)

plt.plot(t, dampedParticle.GetMechanicalEnergy(), label = 'CheckWallLimits clase')
plt.plot(t, dampedParticle2.GetMechanicalEnergy(), label = 'CheckWallLimits modificado')
plt.ylabel('Energía mecánica J(s)')
plt.xlabel('t(s)')
plt.legend()
#plt.show()




print('Tiempo en dejar de rebotar (CheckWallLimits clase): ',  np.round(np.min(t[np.array(dampedParticle.GetMechanicalEnergy()) < -186]), 3), 's' )
print('A 30s, la altura de la particula con (CheckWallLimits modificado)  es ', np.round((dampedParticle2.rVector[-1, 1]+20)/(dampedParticle2.rVector[0, 1] + 20) * 100), r' % de la incial' )

redt = ReduceTime(t, 10, Particles = [dampedParticle])

redt = ReduceTime(t, 10, Particles = [dampedParticle2])

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1)

def init():
    ax.set_xlim(-Limits[0],Limits[0])
    ax.set_ylim(-Limits[1],Limits[1])
    ax.set_xlabel('x(m)')
    ax.set_ylabel('y(m)')

def Update(i, Particles = [], method = 'modificado'):
    
    plot = ax.clear()
    init()
    ax.set_title(r'$t=%.2f \ seconds$ CheckWallLimits ' %(redt[i]) + method, fontsize=15)
    
    for p in Particles:
        x = p.GetRPositionVector()[i,0]
        y = p.GetRPositionVector()[i,1]
     
        
        vx = p.GetRVelocityVector()[i,0]
        vy = p.GetRVelocityVector()[i,1]
        
        circle = plt.Circle( (x,y), p.GetR(), color='k', fill=False)
        ax.add_patch(circle)
        ax.arrow(x,y,vx,vy,color='r',head_width=0.5)
        
    return plot

Animation = anim.FuncAnimation(fig, lambda i: Update(i, Particles = [dampedParticle2]), frames=len(redt),init_func=init, interval = 10, repeat = False)

#plt.show()
plt.close()

"""
