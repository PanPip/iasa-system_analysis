# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
#import pandas as pd
#import collections

#Entering expert data
def calc():
    p = ([[0.65,0.5,0.2,0.5],
          [0.2,0.8,0.65,0.5],
          [0.5,0.35,0.5,0.5],
          [0.8,0.35],
          [0.35,0.65,0.05]])
    
    #Norming the expert scores
    for line in p:
        s = sum(line)
        # #print(s)
        for n, item in enumerate(line):
            line[n] /= s
            
    #to print
    L = ['F1', 'F2', 'F3','F4','F5']
    #p1 = list(p)
    #p1[3].extend((None,None))
    #p1[4].append(None)
    #print(p1)
    '''
    ptab = collections.OrderedDict([ ('Параметр', ['1', '2', '3','4']),
                         ('F1', p1[0]),
                         ('F2', p1[1]),
                         ('F3', p1[2]),
                         ('F4', p1[3]),
                         ('F5', p1[4]) ] )
    df = pd.DataFrame.from_dict(ptab)
    #print(df)    
    '''
    dim = [4,4,4,2,3]
            
    #Entering morf data. Will compose the needed matrixes
    F12 = np.array([[0.5,0.2,0,0.2],
                     [-0.6,0,0,0],
                     [0.5,0.2,0,0.2],
                     [0.3,0,0,0]])
    
    F13 = np.array([[0,0,0,0],
                     [0.7,0,0,0.2],
                     [0,0.5,0.5,0.2],
                     [0,0.3,0.5,0]])
    
    F14 = np.array([[0.7,0.5,0.3,0.5],
                    [-0.3,0.5,0.1,0.3]])
    
    F15 = np.array([[0,0.3,0.5,0.3],
                    [0,0,0,0],
                    [0,-0.2,-0.5,-0.3]])
    
    F23 = np.array([[-0.5,-0.7,-0.3,-0.1],
                     [0.6,0,0.5,0.3],
                     [0,0.7,0,0],
                     [0,0.7,0,0]])
    
    F24 = np.array([[0.4,0,0.5,0.5],
                    [0.4,0.6,0,0]])
    
    F25 = np.array([[-0.2,0.7,-0.3,-0.4],
                    [0,0,0,0],
                    [0.2,-0.7,0.3,0.4]])
    
    F34 = np.array([[-0.1,0.5,0,0],
                    [-0.3,0,0.5,0.7]])
    
    F35 = np.array([[-0.8,-0.4,0.3,0.5],
                    [0,0,0,0],
                    [0.8,0.4,-0.3,-0.5]])
    
    F45 = np.array([[0,0.3],
                    [0,0],
                    [0.3,0]])
    
        
    import pandas
    from functools import reduce
    
    def calc_c(row):
        ar = np.array(row)
        ar = ar[ar>0]
        if len(ar) == 1:
            return float(ar)
        else:
            forw = [2/(1-c)-1 for c in ar]
            mult=1
            for el in forw:
                mult *= el
            res = 1-2/(mult+1)
            return res
        
    
    r = lambda p: np.cos( (np.arccos(1 - 2*p) + np.pi )/3  ) + 1 / 2
    
    def calc_nu(r, p):
        if p >= 0.5:
            res = -np.log2(r(p))
            return res
        else:
            res = - (np.log2(1-r(p))) ** (-1)
            return res
        
    def calc_t(c,nu):
        if nu >= 1:
            res = 1 - 2* ((1-c)/2) ** (nu)
            return res
        else:
            res = 2 * ((1+c)/2) ** (1/nu) - 1
            return res
        
    def calc_p(t):
        res = 3 *((t+1)/2) **2 - 2 * ((t+1)/2) **3
        return res
    
    
    #Function takes the F matr as M
    #also p - the initial probability vector
    #and dim - the dimenstios of the functionl possibilities
    def calc(M,p,dim):
        dim = dim.copy()
        res = []
        row = []
        #Iterating through all possibilities
        for i1 in range(dim[0]):
            for i2 in range(dim[1]):
                for i3 in range(dim[2]):
                    for i4 in range(dim[3]):
                        row.extend((i1+1,i2+1,i3+1,i4+1))
                        for j in range(dim[4]):
                            c_vect = [M[i1][j],M[dim[0]+i2][j],M[dim[1]+i3][j],M[dim[2]+i4][j]] 
                            if not any(c_vect):
                                #Here we go without anything, just return the initial possibility
                                row.append(p[j])
                            else:
                                c = calc_c(c_vect)
                                nu = calc_nu(r,p[j])
                                t = calc_t(c,nu)
                                row.append(calc_p(t))
                        #print(row)
                        res.append(row)
                        row = []
        res = np.array(res)
        return (res)
    
    
    def t_calc(M,p,dim):
        res = []
        row = []
        #Iterating through all possibilities
        for i1 in range(dim[0]):
            for i2 in range(dim[1]):
                row.extend((i1+1,i2+1))
                for j in range(dim[2]):
                    c_vect = [M[i1][j],M[dim[0]+i2][j]] 
                    if not any(c_vect):
                        #Here we go without anything, just return the initial possibility
                        row.append(p[j])
                    else:
                        c = calc_c(c_vect)
                        nu = calc_nu(r,p[j])
                        t = calc_t(c,nu)
                        row.append(calc_p(t))
                        # if i1 == i2== 1:
                            #print (c, p[j],nu,t,calc_p(t))
                #print(row)
                res.append(row)
                row = []
                #print(row)
        res = np.array(res)
        return (res)
    
    M =   np.array([[0,0.5,0.5,0],
                    [0,0,0,0],
                    [0,0.3,0,0],
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,-0.3,0],
                    [0,0,0.3,0],
                    [0.3,0.3,0.3,0],
                    [0.3,0.3,0.3,0],
                    [0.3,0.3,0.3,0],
                    [0,0,0,0.7],
                    [0,0,0,0]])
    
    ptest= ([[0.323,0.516,0.129,0.032],
          [0.342,0.184,0.263,0.053,0.026,0.105,0.026],
          [0.227,0.091,0.364,0.295,0.023]])
    
    
    #res = t_calc (M,ptest[0],[7,5,4])
    #print(res)
    
    def normat(M):
        for row in M:
            s = sum(row[4:])
            for el in range(len(row)-4):
                row[4+el] = row[4+el]/s
                
        return M
    
    def normatest(M):
        for row in M:
            s = sum(row[2:])
            for el in range(len(row)-2):
                row[2+el] = row[2+el]/s
        return M
    
    
    dim = [[4,4,2,3,4],
          [4,4,2,3,4],
          [4,4,2,3,4],
          [4,4,4,3,2],
          [4,4,4,2,3]]
    
    
    #F1
    F1 = np.concatenate((F12,F13,F14,F15), axis=0)
    rF1 = calc (F1,p[0],dim[0])
    nF1 = normat(rF1)
    #F2
    F2 = np.concatenate((np.transpose(F12),F23,F24,F25), axis=0)
    rF2 = calc (F2,p[1],dim[1])
    nF2 = normat(rF2)
    #F3
    F3 = np.concatenate((np.transpose(F13),np.transpose(F23),F34,F35), axis=0)
    rF3 = calc (F3,p[2],dim[2])
    nF3 = normat(rF3)
    #F4
    F4 = np.concatenate((np.transpose(F14),np.transpose(F24),np.transpose(F34),F45), axis=0)
    rF4 = calc (F4,p[3],dim[3])
    nF4 = normat(rF4)
    #F5
    F5 = np.concatenate((np.transpose(F15),np.transpose(F25),np.transpose(F35),np.transpose(F45)), axis=0)
    rF5 = calc (F5,p[4],dim[4])
    nF5 = normat(rF5)
    
    def itermethod(M,p,num):
        p = p.copy()
        M = M.copy()
        pnum = p[num]
        res = np.zeros(len(pnum))
        #print(res)
        del p[num]
        for row in M:
            mult = 1
            for i in range(len(p)):
                #print(row)
                row[i] = p[i][int(row[i])-1]
                mult *= row[i]
            for j in range(len(pnum)):
                #print(mult * row[len(p)+j-1])
                res[j] += mult * row[len(p)+j-1]
        #print(M)
        res[-1] = 0.98 - sum(res[:-1])
        return res    
            
        
    def iterall(M1,M2,M3,M4,M5,p,itr):
        p= p.copy()
        M = [M1,M2,M3,M4,M5]
        for i in range(itr):
            p[0] = list(itermethod(M1,p,0))
            p[1] = list(itermethod(M2,p,1))
            p[2] = list(itermethod(M3,p,2))
            p[3] = list(itermethod(M4,p,3))
            p[4] = list(itermethod(M5,p,4))
            #print(p)
        return p
        
    
    aF1 = itermethod(nF1,p,0)
    aF2 = itermethod(nF2,p,1)
    aF3 = itermethod(nF3,p,2)
    aF4 = itermethod(nF4,p,3)
    aF5 = itermethod(nF5,p,4)
    #print(aF1)
    #print(aF2)
    #print(aF3)
    #print(aF4)
    #print(aF5)
    
    ans = iterall(nF1,nF2,nF3,nF4,nF5,p,100)
    
    #print(ans)
    a1 = ans.copy()
    a1[3].extend((None,None))
    a1[4].append(None)
    '''
    ptab = collections.OrderedDict([ ('Параметр', ['1', '2', '3','4']),
                         ('F1', a1[0]),
                         ('F2', a1[1]),
                         ('F3', a1[2]),
                         ('F4', a1[3]),
                         ('F5', a1[4]) ] )
    '''
    #df = pd.DataFrame.from_dict(ptab)
    #print(df)
    #print(sum(ans[1]))
    #df
    
    def calc_a(M,p,dim):
        dim = dim.copy()
        res = []
        row = []
        #Iterating through all possibilities
        for i1 in range(dim[0]):
            for i2 in range(dim[1]):
                for i3 in range(dim[2]):
                    for i4 in range(dim[3]):
                        for i5 in range(dim[4]):
                            row.extend((i1+1,i2+1,i3+1,i4+1,i5+1))
                            for j in range(dim[5]):
                                c_vect = [M[i1][j],M[dim[0]+i2][j],M[dim[1]+i3][j],M[dim[2]+i4][j],M[dim[3]+i5][j]] 
                                if not any(c_vect):
                                    #Here we go without anything, just return the initial possibility
                                    row.append(p[j])
                                else:
                                    c = calc_c(c_vect)
                                    nu = calc_nu(r,p[j])
                                    t = calc_t(c,nu)
                                    row.append(calc_p(t))
                            #print(row)
                            res.append(row)
                            row = []
        res = np.array(res)
        return (res)
    
    def normat_a(M):
        for row in M:
            s = sum(row[5:])
            for el in range(len(row)-5):
                row[5+el] = row[5+el]/s
                
        return M
    
    def itermethod_a(M,p,num):
        p = p.copy()
        M = M.copy()
        pnum = p[num]
        res = np.zeros(len(pnum))
        #print(res)
        del p[num]
        for row in M:
            mult = 1
            for i in range(len(p)):
                row[i] = p[i][int(row[i])-1]
                mult *= row[i]
            for j in range(len(pnum)):
                #print(mult * row[len(p)+j-1])
                #print(j)
                #print(len(p))
                res[j] += mult * row[len(p)+j-1]
        #print(M)
        return res  
    
    A =    np.array([[0,0.8,0.2,0.4,0.1,0.1],
                     [0.3,0.5,0.7,0.9,0,0.1],
                     [0.2,0.1,0.5,-0.2,0.8,0.6],
                     [0.9,0.1,0.1,-0.2,0.3,0.1],
                     [0.1,0.3,0.6,0.7,0.8,0.3],
                     [0.5,0.6,0.4,0.2,-0.3,0.5],
                     [0,0.4,0.8,0.8,0.6,0.8],
                     [0.2,0.7,0.4,0.3,0,0.7],
                     [-0.8,-0.2,-0.8,-0.3,-0.2,-0.4],
                     [0,0.1,-0.1,0.3,0.3,0.2],
                     [0.3,0.4,0.5,0.3,-0.4,0.4],
                     [0,0.7,0.6,0.4,-0.2,0.2],
                     [0.6,0.7,0.5,0.6,0.8,0.7],
                     [0.2,0.5,0.5,0.7,0,0.3],
                     [-0.3,0,-0.2,-0.1,-0.7,-0.7],
                     [0,0.2,0,0.1,-0.2,0.1],
                     [0.4,0.3,0.4,0.3,0.5,0.5]])
    
    pa = [0.16,0.17,0.17,0.17,0.17,0.16]
    
    dima = [4,4,4,2,3,6]
    
    rA = calc_a (A,pa,dima)
    nA = normat_a (rA)
    np.set_printoptions(threshold=np.nan)
    np.set_printoptions(precision=3)
    #print(nA)
    
    p_a = ans.copy()
    
    p_a.append(pa)
    #print(p_a)
    
    result = list(itermethod_a(nA,p_a,5))
    #print(result)
    return(result)
    
