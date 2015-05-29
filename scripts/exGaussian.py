#!/usr/bin/env python
"""

Exponentially Modified Gaussian
Reference:
Olivier, Jake, and Melissa M. Norberg. "Positively Skewed Data: Revisiting the Box-Cox Power Transformation." International Journal of Psychological Research 3.1 (2010): 68-77.
Note: lamb here corresponds to 1/lambda in the exponential distribution

"""

import numpy as np
import scipy.stats

def pdf(x, params):
    """
    a different formula according to wikipedia 
    (lambda is consistent with the original definition of lambda,
    so have to convert it back before using it
    to keep the consistency of the interface of this module)
    http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    # import scipy.special
    # lamb = 1.0/lamb
    # return lamb/2*np.exp(lamb/2*(2*mu + lamb*sigma*sigma - 2*x)) * scipy.special.erfc((mu + lamb*sigma*sigma - x)/(np.sqrt(2)*sigma))
    """
    mu, sigma, lamb = params
    x = np.array(x)
    return scipy.stats.norm.cdf((x-mu)/sigma - sigma/lamb)*np.exp(sigma*sigma/(2*lamb*lamb)-(x-mu)/lamb)/lamb

def cdf(x, mu, sigma, lamb):
    """
    Reference:
    http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    here lambda is actually lambda, to be consistent with this module, 
    convert lambda back to the normal definition within this function
    """
    lamb = 1.0/lamb
    x = np.array(x)
    u = lamb * ( x - mu )
    v = lamb * sigma
    return scipy.stats.norm.cdf(u,0,v) - np.exp(-u + v*v/2 + np.log(scipy.stats.norm.cdf(u, v*v, v)))
    
def negloglikelihood(*args):
    pdf_vec = pdf(*args)
    log_vec = np.log(pdf_vec[pdf_vec!=0])
    if len(log_vec)==0: return np.inf
    return -np.sum(log_vec)

def mle(x):
    import scipy.optimize
    mu = np.median(x)
    sigma = np.sqrt(np.std(x))
    lamb = 1/np.mean(x)
    params, nll, d = scipy.optimize.fmin_l_bfgs_b(negloglikelihood, [mu, sigma, lamb], fprime=None, args=(x,), approx_grad=True, bounds=[(1e-6, np.inf),(1e-6, np.inf),(1e-6, np.inf)])
    mu, sigma, lamb = params
    print "estimated parames: mu: {0:.2f}, sigma: {1:.2f}, lambda: {2:.2f}".format(mu,sigma, lamb)
    print "negative log likelihood: {0:.0f}".format(nll)
    if d['warnflag'] != 0:
        print "optimization not converged!"
    return mu, sigma, lamb

    
def fit(x):
    n = len(x)
    x = np.array(x)
    m = np.mean(x)
    s = np.std(x)
    g = np.sum(np.power(x - m, 3))/((n-1)*np.power(s,3))
    t = np.power(g/2, 1/3.0)
    mu = m - s * t
    st = 1 - t*t
    if st<0:
        print "skewness greater than 1!"
    #     mu, sigma, lamb = mle(x)
    #     return mu, sigma, lamb
    # else:
    # sigma = s * np.real(np.lib.scimath.sqrt( 1 - t*t ))
    sigma = s * np.sqrt( np.abs(1-t*t))
    lamb = s * t
    return mu,sigma,lamb

def main():
    import matplotlib.pyplot as plt
    print "true params: mu: 10, sigma: 2, lambda: 5"
    x = scipy.stats.norm.rvs(10, 2, 10000)
    y = scipy.stats.expon.rvs(0, 5, 10000)
    z = x + y
    ns,bins,patches = plt.hist(x,50,alpha=0.2, normed=True)
    ns,bins,patches = plt.hist(y,50, alpha=0.2, normed=True)
    ns,bins,patches = plt.hist(z,50, alpha=0.2, normed=True)
    mu, sigma, lamb = fit(z)
    print "esitmated params: mu: {0:.0f}, sigma: {1:.0f}, lambda: {2:.0f}".format(mu,sigma,lamb)
    xp = np.linspace(min(z), max(z), 1000)
    yp = pdf(xp, [mu,sigma,lamb])
    plt.plot(xp, yp, 'm-')
    plt.xlim((0,50))
    plt.show()
    print negloglikelihood(z, (mu,sigma,lamb))
    print negloglikelihood(z, (10,2,5))
    cdf_so_far = yp[0]
    cdf_from_pdf = [ cdf_so_far ]
    for i in xrange(1,len(yp)):
        cdf_so_far += yp[i]
        cdf_from_pdf.append(cdf_so_far)
    cdf_from_pdf = np.array(cdf_from_pdf)/float(cdf_so_far)
    cdf_from_func = cdf(xp, 10, 2, 5)
    plt.plot(xp, cdf_from_pdf, 'b-', lw=1, alpha=0.5)
    plt.plot(xp, cdf_from_func, 'r-', lw=1, alpha = 0.5)
    plt.show()
        
if __name__ == "__main__": main()



