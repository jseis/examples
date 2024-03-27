'''
Created on Mar. 2, 2024

@author: jonah
'''

def effective(r, m):
    #r is original rate
    #m >= 1 is number of compounding terms per year
    #e.g., if r is quoted annually and an equivalent quarterly compounded rate is desired, m = 4
    return (1 + r/m)**m - 1

def accrue_one_per(r, m):
    #discount factor for one period
    #r is annual rate, m is compounding periods per year
    return (1 + r/m)

def accrue_one_year(r, m):
    #discount factor for one year, where rate r is quoted on annual basis and accrue is compounded m times per year
    return accrue_one_per(r, m)**m

def accrue(r, n, m =1):
    return accrue_one_year(r, m)**n

