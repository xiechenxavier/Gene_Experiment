#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 06:21:55 2021

@author: xiechen
"""

import time
import numpy as np
import localsolver
import sys


wij_table = (np.loadtxt("./Gene_Experiment/Gurobi/instances/Gene50_40.txt",delimiter=',')).tolist()

nbGene = len(wij_table)
start = [x for x in range(0, nbGene)] # 0 -> 9
with localsolver.LocalSolver() as ls:
 #
    # Declares the optimization model
    model = ls.model
    # Decision Variable x[i][j]
    x =[[model.bool() for j in range(0,len(start))] for i in range(0,len(start))] 
    n =[model.int(1,len(start)) for i in range(0,len(start))]
    y =[[model.bool() for j in range(0,len(start))] for i in range(0,len(start))] 
    N = len(start)
    # constraints 
    for i in range(0,len(start)):
        for j in range(0,len(start)):
            if(i!=j):
                model.add_constraint(model.eq(x[i][j] + x[j][i],1))
    
    for i in range(0,len(start)):
        model.add_constraint(model.leq(x[i][i],0))
    
    # # constraint 4        
    # for i in range(0,len(start)):
    #     for j in range(0,len(start)):
    #         if(i!=j):
    #             model.add_constraint(model.geq(n[j]+N*(1-x[i][j]),n[i]+1))
        
    # constraint 3
    # for i in range(0,len(start)):
    #     for j in range(0,len(start)):
    #         for k in range(0,len(start)):
    #             if(i!=j!=k):
    #                 model.add_constraint(model.geq(x[i][k], x[i][j]+x[j][k]-1))
     # constraint new4 y[i][j] ==1 => j-i ==1
    for i in range(0,len(start)):
        for j in range(0,len(start)):
            if(i!=j):
                model.add_constraint(model.leq(y[i][j],x[i][j]))
    # constraint new5 任意 (i,j,k) y[i][j] + x[k][i] - x[k][j] <= 1
    for i in range(0,len(start)):
        for j in range(0,len(start)):
            for k in range(0,len(start)):
                if(i!=j!=k):
                    model.add_constraint(model.leq(x[i][j] + y[j][i] + x[j][k]+ x[k][i],2))
   
     # GP3 constraint new6 : y[i][j] + y[j][i] + x[k][i] - x[k][j] <= 1
    # for i in range(0,len(start)):
    #     for j in range(0,len(start)):
    #         for k in range(0,len(start)):
    #             if(i!=j!=k):
    #                 model.add_constraint(model.leq(y[i][j] + y[j][i] + x[k][i] - x[k][j],1))
    
    # # GP1 y[i][j] + x[k][i] - x[k][j] <= 1
    # # GP3 constraint new7 y[i][j] + y[k][j] + y[i][k] + x[k][i] - x[k][j] <=1
    # for i in range(0,len(start)):
    #     for j in range(0,len(start)):
    #         for k in range(0,len(start)):
    #             if(i!=j!=k):
    #                 model.add_constraint(model.leq(y[i][j] + y[k][j] + y[i][k] + x[k][i] - x[k][j],1))
   
    # Maximize value

    MaxProbability = model.sum(wij_table[i][j]*x[i][j] for j in range(0,len(start)) for i in range(0,len(start)))
    
    model.maximize(MaxProbability)
    

    model.close()
    
    # Parameterizes the solver
    if len(sys.argv) >= 4: ls.param.time_limit = int(sys.argv[3])
    else: ls.param.time_limit = 2000


    startslot = time.process_time()
    ls.solve()
    endslot = time.process_time()
    print("it took ",(endslot-startslot)/10,"s to solve")
    
    #
    # treat all solutions such as xij = 1
    # obtenir le classement final
    relation = [[0 for i in range(0,len(start))] for j in range(0,len(start))]
    print(MaxProbability.value)
    for i in range(0,len(start)):
        for j in range(0,len(start)):
            if(i!=j):
                if x[i][j].value:
                    relation[j][i] = 1
                    
    # print(relation)
    sommtab = [sum(l) for l in relation]
    print(sommtab)
    classement = [sommtab.index(i)+1 for i in range(len(start))]
    print(classement)
    
    

                
                



    

       
        
