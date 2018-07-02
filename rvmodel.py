'''
RV models computed with different fitted and fixed parameters.
''' 
from imports import *
import rvcurve

# The standard: fit all parameters
def get_rv1(theta, bjd):
    P,T0,V,K,h,k=theta
    return rvcurve.RV(bjd, theta)

# Fix V
def get_rv2(theta, bjd):
    P,T0,K,h,k=theta
    theta2=P,T0,0.,K,h,k
    return rvcurve.RV(bjd, theta2) 

# Fix period and assume circular orbit
def get_rv3(theta, bjd, P=1.6289):
    T0,V,K=theta
    theta2=P,T0,V,K,0.,0.
    return rvcurve.RV(bjd, theta2)

# Assume everything from Berta15 and T0 from 
# my transit fit 
def get_rv4(theta, bjd, P=1.6289, T0=2457184.55786):
    V,K=theta
    theta2=P,T0,V,K,0.,0.
    return rvcurve.RV(bjd, theta2)

# Assume everything from Berta15 and T0 from 
# my transit fit and no systemic V 
def get_rv5(theta, bjd, P=1.6289, T0=2457184.55786):
    K=theta
    theta2=P,T0,0.,K,0.,0.
    return rvcurve.RV(bjd, theta2)

# Solve for P,T0,K
# assume circular
def get_rv6(theta, bjd):
    P,T0,K=theta
    theta2=P,T0,0.,K,0.,0.
    return rvcurve.RV(bjd, theta2)

# similar to get_rv6 but for 2 planets
def get_rv6_2(theta, bjd):
    P1,T01,K1,P2,T02,K2=theta
    theta1=P1,T01,0.,K1,0.,0.
    theta2=P2,T02,0.,K2,0.,0.
    return rvcurve.RV(bjd, theta1) + rvcurve.RV(bjd, theta2)

# similar to get_rv6_2 but for 3 planets
def get_rv6_3(theta, bjd):
    P1,T01,K1,P2,T02,K2,P3,T03,K3=theta
    theta1=P1,T01,0.,K1,0.,0.
    theta2=P2,T02,0.,K2,0.,0.
    theta3=P3,T03,0.,K3,0.,0.
    return rvcurve.RV(bjd, theta1) + rvcurve.RV(bjd, theta2) + rvcurve.RV(bjd, theta3)
