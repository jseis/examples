'''
Created on Mar. 2, 2024

@author: jonah
'''

from scipy.optimize import fsolve
from metric.interest import accrue_one_per, accrue_one_year

def pres_val(c, r, m = 1):
    #present value
    #c is array of annual cash flows, r is annual rate, m is compounding periods per year
    pv = 0
    discount = 1
    for ci in c:
        pv += ci/discount
        discount *= accrue_one_year(r, m)
    return pv

def ytm(p, c):
    #yield-to-maturity
    #p is current cash flow price
    #c is cash flow per PERIOD. rate is to be given at same period
    #e.g., if cash flow represents bond, c = (c0, c0,..., c0 + F), where F is face value
    def pres_val_diff(r):
        return pres_val(c, r) - p
    return fsolve(pres_val_diff)

def duration(c, r, m = 1):
    #Macauley duration
    #present value times its derivative with respect to r, having accrued one period of interest
    D = 0
    P = pres_val(c, r, m)
    discount = accrue_one_per(r, m)
    for i in range(1, len(c)):
        D += (i/m)*c[i]/discount
        discount *= accrue_one_per(r, m)
    D /= P
    return D

def convexity(c, r, m = 1):
    #present value times its second derivative with respect to r
    C = 0
    P = pres_val(c, r, m)
    discount = accrue_one_per(r, m)
    for i in range(1, len(c)):
        C += i*(i+1)*c[i]/(m**2*discount)
        discount *= accrue_one_per(r, m)
    C /= (P*accrue_one_per(r, m)**2)
    return C
    
def portfolio_pres_val(P):
    #given array P of security present values, calculates present value of portfolio
    return sum(P)

def portfolio_duration(P, D):
    #given array P of security present values and 
    #array D of security durations
    #calculates duration of portfolio
    return (P @ D)/portfolio_pres_val(P)

def portfolio_convexity(P, C):
    #given array P of security present values and 
    #array C of security convexities
    #calculates convexity of portfolio
    return (P @ C)/portfolio_pres_val(P)

def bond_cash_flow(coupon, years, face_value = 1, m = 1):
    #computes a cash flow for a bond given its coupon, face value, 
    #number of payment periods remaining, number of compounding periods per year m
    #default of 1 for face_value corresponds to per-dollar-invested cash flow
    #face_value of x is identifcal to calling with face_value of 1, 
    #then scaling output vector by x
    #c0 is 0
    c = [0]*(years+1)
    for i in range(1,years+1):
        c[i] = coupon/m
    c[years] += face_value
    return c
    

    
    
    