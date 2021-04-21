#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 03:59:14 2021

@author: xiechen
"""
import gurobipy
import numpy as np
import pandas as pd
import time
#导入gurobi库

class Gurobi_Gene:
    
    def __init__(self):
        self.m = gurobipy.Model("GeneSolverModel")
        self.wij_table = []
        
    
    def importInstance(self,filename):
        '''
        cette fonction permet d'importer des donnes d'une fichier
        '''
        #la table de la probabilite d'i est avant d'j
        self.wij_table = (np.loadtxt(filename,delimiter=',')).tolist()
        #la taille est le nombre du Gene
        self.taille = len(self.wij_table)
        self.index=[i for i in range(0,self.taille)]
        self.columns=[j for j in range(0,self.taille)]
        self.packed_wij = pd.DataFrame(self.wij_table,index = self.index, columns = self.columns)
        
        
    def parametrageModel(self):
        '''
        parametrer les contenues communes d'une configuration d'un model
        '''
        #les indices dans un ensemble
        # setPosition = set([ele for ele in range(0,self.taille)])
        # (或者addVars()，如果您希望一次添加多个变量)。变量总是与特定的模型相关联。
        self.x = self.m.addVars(self.index, self.columns, vtype=gurobipy.GRB.BINARY)
        # objectif function
        self.m.setObjective(gurobipy.quicksum(self.x[i, j] * self.packed_wij.at[i, j] 
                                              for i in self.index for j in self.columns),sense=gurobipy.GRB.MAXIMIZE)
        
    # les contraintes principales contiennent constraint 1,2
    def constraintPrincipal(self):
        x = self.x
        for i in range(0,len(self.index)):
            for j in range(0, len(self.columns)):
                if(i!=j):
                    self.m.addConstr(x[i,j]+x[j,i] == 1)
        #constraint 1        
        for i in range(0,len(self.index)):
            self.m.addConstr(x[i,i] <= 0)
            
    # Formulation1: 1,2,4(old)        
    def Formulation1(self):
        #les contraintes principales 1,2
        self.constraintPrincipal()
        #n[i] est la position d'i dans la solution
        self.n = self.m.addVars(self.index,name="n")
        #constraint 4
        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                if(i!=j):
                    self.m.addConstr(self.n[j] + self.taille*(1-self.x[i,j]) >= self.n[i] +1)
                    
    # Formulation2 : 1,2,4,5
    def Formulation2(self):
        #les contraintes principales 1,2
        self.constraintPrincipal()
        
        #varaibles y
        self.y = self.m.addVars(self.index, self.columns, vtype=gurobipy.GRB.BINARY)
        
        #new constraint 4
        for i in range(0,len(self.index)):
            for j in range(0, len(self.columns)):
                if(i!=j):
                    self.m.addConstr(self.y[i,j] <= self.x[i,j])
        #new constraint 5
        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                for k in range(0,len(self.index)):
                    if(i!=j!=k):
                        self.m.addConstr(self.x[i,j] + self.y[j,i] + self.x[j,k]+ self.x[k,i] <= 2)
        
   #Formulation3: 1,2,3,6,7
   # GP : y[i,j] + y[j,i] + x[k,i] - x[k,j] <= 1 (6) y[i,j] + y[k,j] + y[i,k] + x[k,i] - x[k,j] <= 1 (7)
    def Formulation3(self):
        #les contraintes principales 1,2
        self.constraintPrincipal()
        
        #varaibles y
        self.y = self.m.addVars(self.index, self.columns, vtype=gurobipy.GRB.BINARY)
        # positionSet = set(self.index)

        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                for k in range(0,len(self.index)):
                    if(i!=j!=k):
                        # constraint 3
                        self.m.addConstr(self.x[i,j] + self.x[j,k] + 1 - self.x[i,k] <= 2)
                        # new constraint 6
                        self.m.addConstr(self.y[i,j] + self.y[j,i] + self.x[k,i] - self.x[k,j] <= 1)
                        # new constraint 7
                        self.m.addConstr(self.y[i,j] + self.y[k,j] + self.y[i,k] + self.x[k,i] - self.x[k,j] <= 1)
                    
    #Formulation4: 1,2,4,5,6,7 
    def Formulation4(self):
        # Formulation2 : 1,2,4,5
        self.Formulation2()
        
        #ajouter contraintes 6,7
        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                for k in range(0,len(self.index)):
                    if(i!=j!=k):
                        # new constraint 6
                        self.m.addConstr(self.y[i,j] + self.y[j,i] + self.x[k,i] - self.x[k,j] <= 1)
                        # new constraint 7
                        self.m.addConstr(self.y[i,j] + self.y[k,j] + self.y[i,k] + self.x[k,i] - self.x[k,j] <= 1)
                        
    #Formulation5: 1，12，4，5
    def Formulation5(self):
        self.x = self.m.addVars(self.index, self.columns, vtype=gurobipy.GRB.BINARY)
        # new constraint 1
        for i in range(0,len(self.index)):
            for j in range(0, len(self.columns)):
                if(i!=j):
                    self.m.addConstr(self.x[i,j]+self.x[j,i] == 1)
                    
        
        # variables z
        self.z = self.m.addVars(self.index, self.columns,self.index, vtype=gurobipy.GRB.BINARY)
        # variables y
        self.y = self.m.addVars(self.index, self.columns, vtype=gurobipy.GRB.BINARY)
        
        #constraint 15-18
        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                for k in range(0,len(self.index)):            
                    #constraint 15
                    self.m.addConstr(self.z[i,j,k] <= self.x[i,j])
                        
                    #constraint 16
                    self.m.addConstr(self.z[i,j,k] <= self.x[j,k])
                        
                    #constraint 17
                    self.m.addConstr(self.z[i,j,k] <= self.x[i,k])
                        
                    #constraint 18
                    self.m.addConstr(self.z[i,j,k] >= self.x[i,j]+self.x[j,k]-1)
        
        #new constraint 4
        for i in range(0,len(self.index)):
            for j in range(0, len(self.columns)):
                if(i!=j):
                    self.m.addConstr(self.y[i,j] <= self.x[i,j])
        #new constraint 5
        for i in range(0,len(self.index)):
            for j in range(0,len(self.columns)):
                for k in range(0,len(self.index)):
                    if(i!=j!=k):
                        self.m.addConstr(self.x[i,j] + self.y[j,i] + self.x[j,k]+ self.x[k,i] <= 2)
            
     
        
    def PrintResult(self,Formulation):
        #应用哪个求解公式求解
        if(Formulation == 'Formulation1'):
            self.Formulation1()
        elif(Formulation == 'Formulation2'):
            self.Formulation2()
        elif(Formulation == 'Formulation3'):
            self.Formulation3()
        elif(Formulation == 'Formulation4'):
            self.Formulation4()
        elif(Formulation == 'Formulation5'):
            self.Formulation5()
            
        self.m.Params.TimeLimit=2000
        #求解正式开始，测试求解时间
        self.m.optimize()
        
        print("it took ",self.m.Runtime,"s to solve")
         # 输出信息
        dfresult = self.packed_wij * 0
        if self.m.status == gurobipy.GRB.Status.OPTIMAL:
            solution = [k for k, v in self.m.getAttr('x', self.x).items() if v == 1]
            for i, j in solution:
            # print(f"{i} -> {j}：{cost_matrix.at[i,j]}")
                dfresult.at[j, i] = 1
        result = ([[int(ele) for ele in lis]for lis in dfresult.values])
        sommtab = [(sum(l)) for l in result]
        # print(result)
        classement = [sommtab.index(i)+1 for i in range(0,self.taille)]
        return classement
    
    
    def getImportantInfoSolution(self,form):
        if(form == 'Formulation1'):
            self.Formulation1()
        elif(form == 'Formulation2'):
            self.Formulation2()
        elif(form == 'Formulation3'):
            self.Formulation3()
        elif(form == 'Formulation4'):
            self.Formulation4()
        elif(form == 'Formulation5'):
            self.Formulation5()   
        self.m.Params.TimeLimit=300
        #求解正式开始，测试求解时间
        self.m.optimize()
        objective_value = self.m.objVal
        calcul_temps = self.m.Runtime
        return objective_value,calcul_temps
    
    
    def calculAllFormulationTimes(self):
        Formulations = ['Formulation1','Formulation2','Formulation3','Formulation4','Formulation5']
        SolvingInfos = []
        for form in Formulations:
            self.m.reset(0)
            SolvingInfos.append(self.getImportantInfoSolution(form))
        return SolvingInfos
        
if __name__ == '__main__':
    
    Gene_rank1 = Gurobi_Gene()
    Gene_rank1.importInstance("./Gene_Experiment/Gurobi/instances/Gene100_80.txt")
    print(Gene_rank1.wij_table)
    Gene_rank1.parametrageModel()
    print(Gene_rank1.PrintResult("Formulation1"))
    # print(Gene_rank1.calculAllFormulationTimes())
   
    




    
    

