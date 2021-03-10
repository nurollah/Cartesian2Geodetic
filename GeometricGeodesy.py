# In the name of Allah
# -------------------------- Demo module ------------------------------------
# If you have used this module, please also cite the following paper:
# Tatar, N., Farzaneh, S. (2018). Transforming Geocentric Cartesian 
# Coordinates to Geodetic Coordinates by a New Initial Value Calculation Paradigm.
# Journal of the Earth and Space Physics, 44(4), 19-28. 
# doi: 10.22059/jesphys.2018.246251.1006946
# data :  2018
# School of Surveying and Geomatics Engineering, College of Engineering,
# University of Tehran, Iran 
# this code implemented by Nurollah Tatar; Email: n.tatar@ut.ac.ir
# ---------------------------Note ----------------------------------------
# This code is allowed to use only for research purpose, and we
# don't provide any warranty. 
import numpy as np

def Geo2Cart(G,a,e):
    # this code is used to convert geodetic coordinates to cartesian coordinates
    # aboute inputs:
    # G: is a vector that contains [phi, Lamda, Height]
    # a: is the semi-major axis of the biaxial ellipsoid
    # e: is eccentricity of the biaxial ellipsoid
    
    phi =np.deg2rad(G[:,0])
    la = np.deg2rad(G[:,1])
    h = G[:,2]
    m = G.shape[0]
    p=np.zeros([m,3])
    
    N = a/np.sqrt(1-((np.sin(phi))**2)*e**2) # Eq (2).
    p[:,0] = (N+h)*np.cos(phi)*np.cos(la) # X coordinate. Eq (1)
    p[:,1] = (N+h)*np.cos(phi)*np.sin(la) # Y coordinate. Eq (1)
    p[:,2] = (N*(1-e**2)+h)*np.sin(phi) # Z coordinate. Eq (1)

#    for i in range(m):
#        N=a/np.sqrt(1-((np.sin(phi[i]))**2)*e**2) # Eq (2).
#        p[i,0]=(N+h[i])*np.cos(phi[i])*np.cos(la[i]) # X coordinate. Eq (1)
#        p[i,1]=(N+h[i])*np.cos(phi[i])*np.sin(la[i]) # Y coordinate. Eq (1)
#        p[i,2]=(N*(1-e**2)+h[i])*np.sin(phi[i]) # Z coordinate. Eq (1)
    return p

def Cart2Geo(p,a,e):
    # This code is used to convert  cartesian coordinates to geodetic coordinates 
    b = a*np.sqrt(1-e**2)
    e2 = np.sqrt(1-e**2)
    a1 = 1/a
    b1 = 1/b
    a2 = a**2
    b2 = b**2
    
    x = p[:,0]
    y = p[:,1]
    z = p[:,2]
    
    PG = np.sqrt(x**2+y**2) # Eq(5).
    lamda = 2*np.arctan(y/(x+PG)) #  Eq(4).
    lamda = np.rad2deg(lamda)
    lamda = 360*(lamda<0)+lamda
    
    k = np.sqrt((PG*a1)**2+(z*b1)**2) # Eq(9).
    dd = (k-1)*(PG**2+z**2)
    t0 = e2*(k**2*a2+dd)*z/((k**2*b2+dd)*PG+0.000001) # Eq(18).
    C = 1/np.sqrt(e2**2+t0**2) # Eq(20).
    h = (e2*PG+z*t0-b*np.sqrt(1+t0**2))*C # Eq(19).
    N = np.sqrt((a2-(e**2)*(PG-e2*h*C)**2))/e2 # Eq(25).
    phi = np.arctan((N+h)*z/((N*e2**2+h)*PG+0.000001)) # Eq(26).
    phi = np.rad2deg(phi)
    
    G = np.zeros([p.shape[0],3])
    G[:,0] = phi
    G[:,1] = lamda
    G[:,2] = h
    return G
    