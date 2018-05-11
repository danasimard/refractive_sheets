import numpy as np

def __init__():
    return

def mu_and_dtheta(x,thetaR,s,x0):
    dnthetaT = x[0]
    beta = x[1]
    mu0 = x0[0]
    dtheta0 = x0[1]
    
    theta, mus = get_soln(thetaR,dnthetaT,beta,s)
    if (len(theta)==2) and (np.all(theta == theta)):
        dtheta = np.absolute(theta[1]-theta[0])
        mu = np.absolute(mus[0]) #/mu[1] - add back in if we can't assume that the magnification of the los image is 1
    else:
        dtheta = 0.
        mu = 0
        
    return [np.absolute(mu-mu0),np.absolute(dtheta-dtheta0)]

def get_dtheta(beta,thetaR,dnthetaT,s,dtheta_1):
    theta = get_theta(thetaR,dnthetaT,beta,s)
    if (len(theta)==2) and (np.all(theta == theta)):
        dtheta = np.absolute(theta[1]-theta[0])
    else:
        dtheta = 0.
    return dtheta - dtheta_1

def get_theta(thetaR,dnthetaT,beta,s):
    # All parameters, including beta, are dimensionless
    theta = [None,None]
    p = [4.,
         2.*thetaR-8.*beta,
         4.*beta**2.-4.*beta*thetaR,
         2.*beta**2.*thetaR,
         0,
         0,
         -(s*thetaR*dnthetaT)**2.]
    if dnthetaT<0:
        if beta >0:
            theta=np.real(np.sort(np.roots(p)[(np.roots(p)>0)&(np.imag(np.roots(p))==0)])[:2])
            #theta = np.real(np.sort(np.roots(p)))[(np.roots(p)>0)&(np.imag(np.roots(p))==0)]
        else: 
            theta = np.array([np.NaN,beta])
    elif dnthetaT>0:
        theta = np.real(np.sort(np.roots(p)[np.roots(p)>=0])[-1])
        #theta = np.real(np.sort(np.roots(p)[np.roots(p)>=0]))
        if beta < 0:
            theta = np.append(theta,beta)
        else:
            theta = np.append(theta,np.NaN)
    return theta

def get_mu(thetaR,dnthetaT,beta,s):
    theta = get_theta(thetaR,dnthetaT,beta,s)
    mu = (1. / (1. - (s * thetaR**2. * dnthetaT)/( 8. * theta**4. * (1. + thetaR/(2. * theta))**(3./2.)) 
                + ( s * thetaR * dnthetaT)/(theta**3.*(1. + thetaR/(2. * theta))**(1./2.))))
    for i in range(len(mu)): #This isn't working properly
        if theta[i] == beta:
            mu[i] = 1.
    return np.real(mu)

#def get_mu(thetaR,dnthetaT,beta,s):

def get_soln(thetaR,dnthetaT,beta,s):
    theta = get_theta(thetaR,dnthetaT,beta,s)
    mu = (1. / (1. - (s * thetaR**2. * dnthetaT)/( 8. * theta**4. * (1. + thetaR/(2. * theta))**(3./2.)) 
                + ( s * thetaR * dnthetaT)/(theta**3.*(1. + thetaR/(2. * theta))**(1./2.))))
    for i in range(len(mu)): #This isn't working properly
        if theta[i] == beta:
            mu[i] = 1.
    return theta,mu