'''
Created on Mar. 2, 2024

@author: jonah
'''

from metric import interest
from metric import detflow
import numpy as np

if __name__ == '__main__':
    #immunize a portfolio
    #ex. owe 1 million dollars in 7 years
    #use 3 bonds to immunize portfolio against variations in YTM
    #according to 2nd order Taylor expansion of present value
    ytm = 0.08
    n0 = 7
    c0 = detflow.bond_cash_flow(0, 7) #obligation is same as zero coupon bond
    P0 = 1/interest.accrue(ytm, n0)
    D0 = 0
    C0 = 0
    
    #bond 1
    n1 = 6
    r1 = 0.07
    c1 = detflow.bond_cash_flow(r1, n1) #bond cash flow 1, normalized to face value of 1
    D1 = detflow.duration(c1, ytm)
    C1 = detflow.convexity(c1, ytm)
    
    #bond 2
    n2 = 6
    r2 = 0.1
    c2 = detflow.bond_cash_flow(r2, n2) #bond cash flow 2, normalized to face value of 1
    D2 = detflow.duration(c2, ytm)
    C2 = detflow.convexity(c2, ytm)
    
    #bond 3
    n3 = 9
    r3 = 0.02
    c3 = detflow.bond_cash_flow(r3, n3) #bond cash flow 3, normalized to face value of 1
    D3 = detflow.duration(c3, ytm)
    C3 = detflow.convexity(c3, ytm)
    
    M = np.array([[1, 1, 1], [D1/P0, D2/P0, D3/P0], [C1/P0, C2/P0, C3/P0]])
    b = np.array([P0, D0, C0])
    P = np.linalg.solve(M, b) #amount invested in each bond
    
    print("P1 = " + str(P[0]))
    print("P2 = " + str(P[1]))
    print("P3 = " + str(P[2]))
    
    
    
    
    
    
    
    
    