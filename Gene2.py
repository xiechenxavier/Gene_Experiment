#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 13:47:18 2021

@author: xiechen
"""
import numpy as np
from itertools import product
import localsolver
import sys

wij_table = (np.loadtxt("./Gene/Gurobi/instances/Gene50_40.txt",delimiter=',')).tolist()
nb_genes = len(wij_table)
start = [x for x in range(0, nb_genes)] # 0 -> 9

with localsolver.LocalSolver() as ls:
    
    model = ls.model
    genes = model.list(nb_genes)
    
    model.constraint(model.count(genes) == nb_genes)
    
    
    estAvant = model.array(wij_table)
    
    #maximize
    selector = model.lambda_function(lambda i,j: model.at(estAvant,genes[i],genes[j]) 
                           if  j>i else 1-model.at(estAvant,genes[i],genes[j]))
    #                              
    obj = model.sum((model.range(0,nb_genes),5),selector)
    
    # obj = 0
    # for i in range(0,nb_genes-2):
    #     for j in range(i+1,nb_genes-1): 
    #         obj += model.at(estAvant,genes[i],genes[j])
    
    model.maximize(obj)
    
    model.close()
    
    if len(sys.argv) >= 4: ls.param.time_limit = int(sys.argv[3])
    else: ls.param.time_limit = 200
    
    ls.solve()
    # print(obj.value)
    # for g in genes.value:
    #     print(g)
    print(obj.value)
