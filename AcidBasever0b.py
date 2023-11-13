# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 16:50:03 2019

@author: Bert
"""

import pandas as pd
from numpy import exp
from scipy.optimize import brentq


species_par = pd.read_excel(r'C:\Users\GBia\Source\Repos\Qingdian-Shu\ConsenCUS-BPMED\BasicParSpeciesBook.xlsx', sheet_name='Sheet1',index_col=0)

Fa = 96485
Ts = 298
Rg =8.314

def pK(name_species):
    return species_par.at[name_species,'pK']

def q(name_species):
    return species_par.at[name_species,'charge']

def Dfac(name_species):
    return species_par.at[name_species,'Dfactor']

def Phase(name_species):
    return species_par.at[name_species,'phase']

pKw = pK('H2O')


class Compound():
#species is a number of species names that are connect via acid/base reaction, 
#e.g. carbonic acid and carbonates,ordered via pron donating  order 
# gasouos compound should the last in the species list    
# e.g.CO2g, H2CO3, HCO3 CO3    
# a number of methods have been defined
    
    def __init__(self,name_compound,species_all):
        self.name = name_compound
        self.species = species_all
        self.par = species_par.loc[self.species]
        self.n_species = len(self.species)
        if Phase(self.species[-1]) == 'g': 
            self.volatile = 1  
            for i in range(self.n_species-1):
                if q(self.species[i]) == 0:
                    self.vol_species = i
        else: 
            self.volatile = 0
            self.vol_species = 'NaN'
        self.n_solutes =self.n_species-self.volatile
        return
                        
    def fractions(self,pH):
        #calculates the fractions of the different aq species at a given pH
        H = 10**(-pH)
        frac =[1]
        for i in range(self.n_solutes-1):     
            frac.append(frac[i]*10**(-pK(self.species[i]))/H)
        sum_frac = sum(frac)                
        for i in range(self.n_solutes):
            frac[i]= frac[i]/sum_frac
        return  frac       

    def ct_from_p(self,pH,p):
        #calculates the total concetration for open system with given pH and pressure p
        frac = self.fractions(pH)
        c_volatile = 10**(-pK(self.species[-1]))*p
        c_t = c_volatile/frac[self.vol_species] 
        return c_t
        
    def state(self,open_sys,pH,pC):
        #calculates the concetration of all species given pH and
        #if open pC is the gas pressure while if not open it pc is the total amount
        if (open_sys == 1 and self.volatile == 1) :
            c_tot = self.ct_from_p(pH,pC) 
        else: 
            c_tot = pC
        frac = self.fractions(pH)                   
        c = [i for i in range(self.n_solutes)]   #define list for state C
        for i in range(self.n_solutes):        # calculate c now we know c total
                c[i]= frac[i]*c_tot
        return c
        
    def charge_total(self, open_sys,pH,pC):
        #calculate the total charge of the system
        c = self.state(open_sys,pH,pC)
        charge_tot = 0
        for i in range(self.n_solutes):
             charge_tot =charge_tot +q(self.species[i])*c[i]
        return charge_tot
    
    
    def concentration_total(self, open_sys,pH,pC):
        # calculate the total charge of the compound
        c = self.state(open_sys,pH,pC)
        c_tot = 0
        for i in range(self.n_solutes):
                  c_tot = c_tot + c[i]
        return c_tot
    
    def partial_conductivity(self,open_sys,pH,pC,eps):
        #calculates the partial conductivity of compound
        # eps = 10**-5 gives conductivity as mS/cm
        f = (Fa**2)/(Rg*Ts)
        c = self.state(open_sys,pH,pC)
        l = 0
        for i in range(self.n_solutes):
             l = l +(q(self.species[i])**2)*Dfac(self.species[i])*c[i]
        return l*f*eps
    
    def p_from_ct(self,open_sys,pH,pC):
        #calculates the total concetration for open system with given pH and pressure p
        c = self.state(open_sys,pH,pC)
        p= c[self.vol_species]/10**(-pK(self.species[-1]))
        return p
                  
    def show_all(self):
        print(self.species)
        print('-------')
        print(self.par)
        print('-------')
        print('num species = ', self.n_species,'  volatile yes/no =  ',self.volatile, 'volatile species =  ',self.vol_species)
        return

class Solution():
    #compounds is list of compounds
    #condition is list of values of Ct/pG order of structure
    #water is alsway implitely defined with solution and global pKw
    
    def __init__(self,name_solution,compounds_all,balance_values,open_sys =0):
        self.name = name_solution
        self.n_compounds = len(compounds_all)
        self.bal = balance_values
        self.sys_type = open_sys
        self.comp = compounds_all
        return 
    
    def update(self,balance):
        self.bal = balance
        return
    
    def update_type(self,open_sys):
        self.sys_type = open_sys
        return
           
    def charge_balance(self,pH):
        #caluclate the charge balance based on the compounds charge
        tot_charge = 10**(-pH)-10**(pH-pKw)
        for i in range(self.n_compounds):
            tot_charge = tot_charge + self.comp[i].charge_total(self.sys_type,pH,self.bal[i])
        return tot_charge

    def equi(self):
        # find pH at which charge balance = 0
        y = brentq(self.charge_balance,0,15)
        return y       

    def state(self):
        #calculate the concentration of all solutiom species
        pHe = self.equi()
        self.solution_state = []
        for i in range(self.n_compounds):
            i_state = self.comp[i].state(self.sys_type,pHe,self.bal[i])
            for j in range(self.comp[i].n_solutes):
                self.solution_state.append([self.comp[i].species[j] ,i_state[j]])
        self.solution_state.append(["H",10**(-pHe)])
        self.solution_state.append(["OH",10**(pHe-pKw)])
        return self.solution_state
    
    def conductivity(self,eps):
        pHe = self.equi()
        f = (Fa**2)/(Rg*Ts)
        y = 0
        for i in range(self.n_compounds):
            y = y + self.comp[i].partial_conductivity(self.sys_type,pHe,self.bal[i],eps)
        y = y + (q('H')**2)*Dfac('H')*10**(-pHe)*eps*f
        y = y + (q('OH')**2)*Dfac('OH')*10**(pHe-pKw)*eps*f
        return y 
          
    def characteristics(self):
        #calculate  some specific numvers for solution ph and total of compunds
        pHe = self.equi()
        self.solution_totals =[]
        self.solution_totals.append(['pH', pHe])   
        for i in range(self.n_compounds):
            if self.comp[i].volatile == 1:
                p_i = self.comp[i].p_from_ct(self.sys_type,pHe,self.bal[i])
                names_i = self.comp[i].species[-1]
                self.solution_totals.append([names_i, p_i])
        for i in range(self.n_compounds):
            names_i = self.comp[i].name
            c_i = self.comp[i].concentration_total(self.sys_type,pHe,self.bal[i])
            self.solution_totals.append([names_i, c_i])
        self.solution_totals.append(['Conductivity', self.conductivity(10**-5)])    
        self.solution_totals.append(['EN', self.charge_balance(pHe)])
        return self.solution_totals
         
    def save_excell(self,name, show_on = False):
        writer = pd.ExcelWriter(name + '.xlsx', engine='xlsxwriter')
        frame_temp = pd.DataFrame(self.characteristics(),columns=['name','mol/L or [1]'])
        if show_on :print(frame_temp)
        frame_temp.to_excel(writer,index = False)
        frame_temp = pd.DataFrame(self.state(),columns=['name','mol/L or [1]'])
        if show_on :print(frame_temp)
        frame_temp.to_excel(writer,index = False, startcol = 3)
        z = [['Open Sys?',self.sys_type]]
        for i in range(self.n_compounds):
            z.append([self.comp[i].name,self.bal[i]])
        frame_temp = pd.DataFrame(z,columns=['name','mol/L or [1]'])
        if show_on :print(frame_temp)
        frame_temp.to_excel(writer,index = False, startrow = 8)
        writer._save()
        return
# next functions needed to describe the composition of resin in equi with solution
# this will be denoted the associated resin       
    def K_res(self,E):
        E_dim = E*Fa/(Rg*Ts)
        K = exp(E_dim)
        return K
        
    def ass_res_charge_balance(self,E,res_comp,x_fixed):
        tot_charge = q(res_comp.species[0])*x_fixed
#        print(tot_charge,'  ',self.K_res(E))
        n_species = len(self.state())  
        for i in range(n_species):
            tot_charge = tot_charge + (self.K_res(E)**(-q(self.state()[i][0]))*self.state()[i][1]*q(self.state()[i][0]))
#            print(i,' ',tot_charge,'    ',(self.state()[i][0]),'    ',self.state()[i][1],'  ',q(self.state()[i][0]))
        return tot_charge    

    def equi_E_res(self,res_comp,x_fixed):
        y = brentq(lambda x: self.ass_res_charge_balance(x,res_comp,x_fixed),-1,1)    
        return y    
    
    def ass_res_characteristics(self,res_comp,x_fixed):
        E = self.equi_E_res(res_comp,x_fixed)
#        print('E = ',E)
        Kd = self.K_res(E)
        pH_sol = self.equi()
        res_con =[]
        for i in range(self.n_compounds):
            c_t  =0
            fractions = self.comp[i].fractions(pH_sol)
            for j in range(self.comp[i].n_solutes): 
                c_t = c_t + (Kd**(-q(self.comp[i].species[j])))*fractions[j]
            if (self.comp[i].volatile == 1 and self.sys_type ==1):
                c_t = c_t*self.comp[i].ct_from_p(pH_sol,self.bal[i])
            else:
                c_t = c_t*self.bal[i]
            res_con.concat([self.comp[i].name,c_t])            
        res_con.concat([res_comp.name, x_fixed])    
        res_con.concat(['E', E ])      
        res_con.concat(['EN', self.ass_res_charge_balance(E,res_comp,x_fixed)]) 
        return res_con 
    
    def ass_res_compounds(self,res_comp,x_fixed):
        z = self.ass_res_characteristics(res_comp,x_fixed)
        comp_list,bal_list = [],[]
        for i in range(self.n_compounds):
            comp_list.append(self.comp[i])
            bal_list.append(z[i][1])
        comp_list.append(res_comp)
        bal_list.append(x_fixed)
        return comp_list,bal_list    

# next functions needed to describe the composition of solution  in equi with resion solution
# we refer to this as the associated solution
        
    def get_n_resin(self):
        # calculates the location of the resin species in the state vector
        n_resin =-1
        n_species = len(self.state()) 
        for i in range(n_species):
            if Phase(self.state()[i][0])== 'r':
                n_resin =i
        return n_resin
    
    def ass_sol_charge_balance(self,E,n_resin):
        n_species = len(self.state())  
        tot_charge=0
        for i in range(n_species):
            if i != n_resin:
               tot_charge = tot_charge + (self.K_res(E)**(-q(self.state()[i][0]))*self.state()[i][1]*q(self.state()[i][0]))    
        return tot_charge
  
    def equi_E_sol(self):
        n_resin = self.get_n_resin()
        y = brentq(lambda x: self.ass_sol_charge_balance(x,n_resin),-1,1)    
        return y         
  
    def ass_sol_characteristics(self):
        n_resin = self.get_n_resin()
        E = self.equi_E_sol()
#        print('E = ',E)
        Kd = self.K_res(E)
        pH_sol = self.equi()
        res_con =[]
        for i in range(self.n_compounds):
            c_t  =0
            fractions = self.comp[i].fractions(pH_sol)
            for j in range(self.comp[i].n_solutes): 
                c_t = c_t + (Kd**(-q(self.comp[i].species[j])))*fractions[j]
            c_t = c_t*self.bal[i] 
            if Phase(self.comp[i].species[0]) != 'r':
               res_con.concat([self.comp[i].name,c_t])            
        res_con.concat(['E', E ])      
        res_con.concat(['EN', self.ass_sol_charge_balance(E,n_resin)]) 
        return res_con 

    def ass_sol_compounds(self):
        z = self.ass_sol_characteristics()
        comp_list,bal_list = [],[]
        j=0
        for i in range(self.n_compounds):
            if Phase(self.comp[i].species[0]) != 'r':
                bal_list.append(z[j][1])
                comp_list.append(self.comp[j])
                j = j+1
        return comp_list, bal_list


def compound_examples():
    Ct = Compound('Ctot',['H2CO3','HCO3','CO3','CO2g'])
    Nt = Compound('Nt',['NH4','NH3'])
    AnExt = Compound('AnExt',['AnEx'])
    print('Ct.show_all()'.center(80,'+'))
    Ct.show_all()
    print('Nt.show_all'.center(80,'+'))
    Nt.show_all()
    print('AnExt.show_all()'.center(80,'+'))
    AnExt.show_all()
    print('state(open_sys,pH,pC) par = 1,7,1 is thus open'.center(80,'+'))
    print(Ct.state(1,7,1))
    print('state(open_sys,pH,pC) par = 0,7,1, is thus closed'.center(80,'+'))
    print(Ct.state(0,7,1))
    print('state(open_sys,pH,pC) par = 1,6.35 ,1'.center(80,'+'))
    print(Ct.state(1,6.35,1))
    print('Ct.charge_total(open_sys,pH,pC) par1,6.35 ,1'.center(80,'+'))
    print(Ct.charge_total(1,6.35,1))
    print('Ct.concentration_total(open_sys,pH,pC) par1,6.35 ,1'.center(80,'+'))
    print(Ct.concentration_total(1,6.35,1))
    print('Ct.partial_conductivity(open_sys,pH,pC,eps) par1,6.35 ,1'.center(80,'+'))
    print(Ct.partial_conductivity(1,6.35,1,10**-5))
    
def main():       
    Ct = Compound('Ctot',['H2CO3','HCO3','CO3','CO2g'])
    Nt = Compound('Nt',['NH4','NH3'])
    AnExt = Compound('AnExt',['AnEx'])
    Nat = Compound('Na',['Na'])
    sol = Solution('test',[Nt,Nat,Ct],[0.075,0.025,1],1)
    print(sol.name.center(80,'+'))
    print(*sol.state(),sep ='\n')
    print('characteristics'.center(80,'+'))
    print(*sol.characteristics(),sep= '\n')
    sol.save_excell(sol.name)
    print('E resin '.center(80,'+'))
    print(sol.equi_E_res(AnExt,2))
    print('ass_res_charateristics '.center(80,'+'))
    print(*sol.ass_res_characteristics(AnExt,2),sep ='\n')
    y=sol.ass_res_compounds(AnExt,2)
    res= Solution('resin',y[0],y[1],open_sys=0)
    print(res.name.center(80,'+'))
    print(*res.characteristics(),sep = '\n')
    print('state'.center(80,'+'))
    print(*res.state(),sep ='\n')
#    for i in range(len(sol.state())):
#        print(sol.state()[i][0],'  ',sol.state()[i][1]/sol_res.state()[i][1])
    print('resin associated solution characteristics '.center(80,'+'))
    print(*res.ass_sol_characteristics(),sep = '\n')
    z = res.ass_sol_compounds()
    ass_sol = Solution('associated solution',z[0],z[1],open_sys=0)
    print('resin associated solution state '.center(80,'+'))
    print(*ass_sol.state(),sep = '\n')
    return 
 
#compound_examples()
