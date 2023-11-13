# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:50:17 2019

This program describes the performance given a certain transfer efficiency of 
the potassium. 

@author: Bert
"""

from re import A, L
import AcidBasever0b as ab
import pandas as pd
import frameplot as fp
import time
from numpy import log
from numpy import linspace
from numpy import ndarray
from scipy.optimize import brentq

An_cond = pd.read_excel(r'C:\Users\GBia\Source\Repos\Qingdian-Shu\ConsenCUS-BPMED\Conductivity calculation BPMED 1.1M.xlsx', sheet_name='Analyte Conductivity Matrix',index_col='Index')
Cat_cond = pd.read_excel(r'C:\Users\GBia\Source\Repos\Qingdian-Shu\ConsenCUS-BPMED\Conductivity calculation BPMED 1.1M.xlsx', sheet_name='Catholyte Conductivity Matrix',index_col='Index')

class ed_cell_ideal():
    
    def __init__(self,K_in,Ct_in,anions,anions_in,cell_charac):
        self.q_area = cell_charac[0]
        self.l_an = cell_charac[1]
        self.l_m = cell_charac[2]
        self.l_cat = cell_charac[3]
        self.l_tot = cell_charac[4]
        self.mem_c = cell_charac[5] 
        self.pCO2_an = cell_charac[6]
        self.Ct = ab.Compound('Ctot',['H2CO3','HCO3','CO3','CO2g'])
        self.Kt =ab.Compound('K',['K'])
        self.an_in=anions_in
        self.sol_comp=[self.Kt,self.Ct]+anions
        self.sol_1 = ab.Solution('Solution_in',self.sol_comp,[K_in,Ct_in]+self.an_in,0)
        return

    def get_I_tot(self, K_eff):
        # sol_1 is the effluent of the absorber defined in __init__
        # sol_2 is the product of mixing the recycle with the influent
        # sol_3 is the solution after the gas exchange
        # sol_4 is the solution after acidification
        # sol_5 is the catholyte 
        alpha = 4.75
        lamb_r = 100
        
        ########################define all the solutions with internal degasing########################
        s1 = self.sol_1.characteristics()
        K_1 = s1[2][1] 
        K_2 = K_1*(1+lamb_r-lamb_r*K_eff)/(1+lamb_r)
        K_3 = K_1 * (1-K_eff)
        K_4 = K_3
        K_5 = K_1
        Ct_1 =s1[3][1]
        sol_4 = ab.Solution('sol_4',self.sol_comp,[K_4,self.pCO2_an]+self.an_in,1)
        Ct_4 = sol_4.characteristics()[3][1]
        if Ct_4 < Ct_1:
            self.sol_4 = sol_4
        else:
            self.sol_4 = ab.Solution('sol_4',self.sol_comp,[K_4,Ct_1]+self.an_in,0)
        s4 = self.sol_4.characteristics()
        Ct_4 = s4[3][1]
        Ct_2 = (Ct_1+lamb_r*Ct_4)/(1+lamb_r)
        sol_2 = ab.Solution('sol_2',self.sol_comp,[K_2,Ct_2]+self.an_in,0)
        pCO2_2 = sol_2.characteristics()[1][1]
        if pCO2_2 > self.pCO2_an:
            self.sol_2 = ab.Solution('sol_2',self.sol_comp,[K_2,self.pCO2_an]+self.an_in,1)
        else:
            self.sol_2 = sol_2
        s2 = self.sol_2.characteristics()
        Ct_2 = s2[3][1]
        #sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,self.pCO2_an]+self.an_in,1)
        #Ct_3 = sol_3.characteristics()[3][1]
        #if Ct_3 < Ct_2:
        #    self.sol_3 = sol_3
        #else: 
        #    self.sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,Ct_2]+self.an_in,0)
        self.sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,Ct_2]+self.an_in,0)
        s3 = self.sol_3.characteristics()
        Ct_3 = s3[3][1]
        Ct_5 = Ct_4
        self.sol_5 = ab.Solution('sol_5',self.sol_comp,[K_5,Ct_5]+self.an_in,0)
        s5 = self.sol_5.characteristics()
        an = self.sol_3.state()
        cat = self.sol_5.state()
        ##############################################################################################
        
        I_K = self.q_area*K_1*K_eff
        f_Ha = an[-2][1]/(an[0][1]+an[-2][1])
        f_Hc = cat[-2][1]/(cat[0][1]+cat[-2][1])
        f_m = (f_Ha+f_Hc)/2
        d_H = f_Hc-f_Ha
        I_tot = (I_K-ab.Dfac('K')*10**-8*self.mem_c/self.l_m*d_H)*(1+(alpha-1)*f_m)/(1-f_m)-(alpha-1)*ab.Dfac('K')*10**-8*self.mem_c/self.l_m*d_H
        return I_tot

    def get_K_eff(self, I_dens):
        # calculates the K effciency that belongs to a certain current
        I_dif = lambda x:self.get_I_tot(x)-I_dens/ab.Fa
        z= brentq(I_dif,0,0.99999)
        return z
          
    def all_state(self, I_dens):
        # determines the composition of the anode
        alpha = 4.75
        lamb_r = 100

        ########################define all the solutions with external degasing########################
        s1 = self.sol_1.characteristics()
        K_1 = s1[2][1]
        K_eff = self.get_K_eff(I_dens)
        K_2 = K_1*(1+lamb_r-lamb_r*K_eff)/(1+lamb_r)
        K_3 = K_1 * (1-K_eff)
        K_4 = K_3
        K_5 = K_1
        Ct_1 =s1[3][1]
        sol_4 = ab.Solution('sol_4',self.sol_comp,[K_4,self.pCO2_an]+self.an_in,1)
        Ct_4 = sol_4.characteristics()[3][1]
        if Ct_4 < Ct_1:
            self.sol_4 = sol_4
        else:
            self.sol_4 = ab.Solution('sol_4',self.sol_comp,[K_4,Ct_1]+self.an_in,0)
        s4 = self.sol_4.characteristics()
        Ct_4 = s4[3][1]
        Ct_2 = (Ct_1+lamb_r*Ct_4)/(1+lamb_r)
        sol_2 = ab.Solution('sol_2',self.sol_comp,[K_2,Ct_2]+self.an_in,0)
        pCO2_2 = sol_2.characteristics()[1][1]
        if pCO2_2 > self.pCO2_an:
            self.sol_2 = ab.Solution('sol_2',self.sol_comp,[K_2,self.pCO2_an]+self.an_in,1)
        else:
            self.sol_2 = sol_2
        s2 = self.sol_2.characteristics()
        Ct_2 = s2[3][1]
        #sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,self.pCO2_an]+self.an_in,1)
        #Ct_3 = sol_3.characteristics()[3][1]
        #if Ct_3 < Ct_2:
        #    self.sol_3 = sol_3
        #else: 
        #    self.sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,Ct_2]+self.an_in,0)
        self.sol_3 = ab.Solution('sol_3',self.sol_comp,[K_3,Ct_2]+self.an_in,0)
        s3 = self.sol_3.characteristics()
        Ct_3 = s3[3][1]
        Ct_5 = Ct_4
        self.sol_5 = ab.Solution('sol_5',self.sol_comp,[K_5,Ct_5]+self.an_in,0)
        s5 = self.sol_5.characteristics()
        an = self.sol_3.state()
        cat = self.sol_5.state()
        ##############################################################################################

        print(*s1,sep= ' \ ')
        #print(Ct_1,'   ',Ct_2,' ',Ct_3)
        # print('solution 2 charateristics')
        print(*s2,sep ='\n')
        # print('solution 2 state')
        # print(*self.sol_2.state(),sep ='\n')
        # print(*self.sol_3.state(),sep ='\n')
        #I_K = self.q_area*K_1*K_eff
        f_Ha = an[-2][1]/(an[0][1]+an[-2][1])
        f_Hc = cat[-2][1]/(cat[0][1]+cat[-2][1])
        f_m = (f_Ha+f_Hc)/2
        d_H = f_Hc-f_Ha
        #I_tot = (I_K-ab.Dfac('K')*self.mem_c/self.l_m*d_H)*(1+(alpha-1)*f_m)/(1-f_m)-(alpha-1)*ab.Dfac('K')*self.mem_c/self.l_m*d_H
        I_tot = I_dens/ab.Fa
        I_K = (I_tot + ab.Dfac('K') * 10**-8 * self.mem_c / self.l_m * (alpha-1) * d_H) * (1 - f_m) / (1 + (alpha - 1) * f_m) + ab.Dfac('K') * 10**-8 * self.mem_c / self.l_m * d_H
        #K_eff = I_K / self.q_area / K_1
        t_K = I_K/I_tot
        #e_an = self.l_an*I_tot*ab.Fa/(s2[-2][1]*0.01)
        #e_mem = 0.025*log(s3[2][1]/s2[2][1])
        #e_cat = self.l_cat*I_tot*ab.Fa/(s3[-2][1]*0.01)
        #e_reac = 0.059*(s3[0][1]-0.29)
        #power = (e_an+e_mem+e_cat+e_reac)*I_tot*ab.Fa
        #flow_CO2 = self.q_area*(Ct_1-Ct_2)
        #power_spec= power/flow_CO2/1000

        ####################Find the accurate conductivity from OLI simulation######################
        Ct_5_OLI = round(s5[3][1], 3)
        K_3_OLI = round(s3[2][1], 3)
        if K_3_OLI > 2.100:
            K_3_OLI = 2.100
        if K_3_OLI < 0:
            K_3_OLI = 0
        if Ct_5_OLI > 2.000:
            Ct_5_OLI = 2.000
        s3_cond_np = An_cond.at[K_3_OLI, 'Conductivity (mS/cm)']
        s3_cond = s3_cond_np.astype(float) #conductivity of acidifying solution without gas bubbles
        s5_cond_np = Cat_cond.at[Ct_5_OLI, 'Conductivity (mS/cm)']
        s5_cond = s5_cond_np.astype(float) #conductivity of catholyte without gas bubbles
        #calculate gas void fraction
        s3_CO2_gas = s1[3][1] - s4[3][1]
        CO2_void = s3_CO2_gas * ab.Rg * ab.Ts / (s3_CO2_gas * ab.Rg * ab.Ts + self.pCO2_an * 101.325 * (1 + lamb_r)) #derived from Vg/(Vg+Vl)
        if CO2_void <= 0.12:
            s3_cond_corrected = s3_cond * (1 - CO2_void) ** 1.5
        else:
            s3_cond_corrected = s3_cond * (1 - CO2_void) / (1 + CO2_void / 2)
        e_an = self.l_an*I_tot*ab.Fa/(s3_cond_corrected*0.01)
        e_mem = 0.025*log(s5[2][1]/s3[2][1])
        #e_cat = self.l_cat*I_tot*ab.Fa/(s5_cond*0.01)
        #e_reac = 0.059*(s3[0][1]-s2[0][1])
        #e_reac = 0.059*(s5[0][1]-0.29)
        #power = (e_an+e_mem+e_cat+e_reac)*I_tot*ab.Fa
        #flow_CO2 = self.q_area*(Ct_1-Ct_4)
        #power_spec= flow_CO2*1000/power
        #power_spec= power/flow_CO2/1000
        ############################################################################################
        
        ####################Include H2 void fraction######################
        #s5_H2_gas = I_dens / ab.Fa / 2 / (self.q_area * lamb_r)
        #H2_void = s5_H2_gas * ab.Rg * ab.Ts / (s5_H2_gas * ab.Rg * ab.Ts + 0.7 * 101.325) #derived from Vg/(Vg+Vl)
        #if H2_void <= 0.12:
            #s5_cond_corrected = s5_cond * (1 - H2_void) ** 1.5
        #else:
            #s5_cond_corrected = s5_cond * (1 - H2_void) / (1 + H2_void / 2)
        e_cat = self.l_cat*I_tot*ab.Fa/(s5_cond*0.01)
        #e_reac = 0.059*(s3[0][1]-s2[0][1])
        e_reac = 0.059*14
        power = (e_an+e_mem+e_cat+e_reac)*I_tot*ab.Fa
        flow_CO2 = self.q_area*(Ct_1-Ct_4)
        #power_spec= flow_CO2*1000/power
        power_spec= power/flow_CO2/1000
        ##################################################################

        #creating res 
        res = []
        res.append(K_eff)
        res.append(t_K)
        res.append(d_H)
        #res.append(K_eff/t_K) #*(K_1/(K_1-s1[4][1])))
        res.append(I_dens/K_1/self.q_area/ab.Fa)
        res.append(1-(Ct_4/Ct_1))
        res.append((Ct_1-Ct_4)/K_1*t_K/K_eff)
        res.append(I_tot*ab.Fa)
        res.append(flow_CO2)
        res.append(e_an)
        res.append(e_mem)
        res.append(e_cat)
        res.append(e_reac)
        res.append(power)
        res.append(power_spec)
        res.append(s1[3][1]) # 3= Ct
        res.append(s2[3][1])
        res.append(s3[3][1])
        res.append(s4[3][1])
        res.append(s5[3][1])
        res.append(s1[0][1]) # 0 =pH
        res.append(s2[0][1])
        res.append(s3[0][1])
        res.append(s4[0][1])
        res.append(s5[0][1])
        #res.append(s1[-2][1]) # -2 is conductivity
        #res.append(s2[-2][1])
        #res.append(s3[-2][1])
        res.append(s3_cond)
        res.append(s5_cond)
        res.append(s1[1][1]) # 1 is pCO2g
        res.append(s2[1][1])
        res.append(s3[1][1])
        res.append(s4[1][1])
        res.append(s5[1][1])
        res.append(s1[2][1]) # 2 is K
        res.append(s2[2][1])
        res.append(s3[2][1])
        res.append(s4[2][1])
        res.append(s5[2][1])
        res.append(CO2_void)
        res.append(s3_cond_corrected)
        return res
    
    def trace(self, L_r_min = 0.1, L_r_max = 1.2, num_step = 100, I_dens = 1):
        if L_r_min <= 0:
            L_r_min = 0.1
        load_range = linspace(L_r_min, L_r_max, num_step)
        s1 = self.sol_1.characteristics()
        K_1= s1[2][1] 
        self.q_area = I_dens / (L_r_min * K_1 * ab.Fa)
        res = self.all_state(I_dens) 
        cols=['K_eff','t_K','fH','Load_R','Ct_rem','Ct_coul_eff','I_dens','J_CO2','E_an','E_m','E_cat','E_reac','power','Co2 power','Ct_1','Ct_2','Ct_3','Ct_4','Ct_5','pH_1','pH_2','pH_3','pH_4','pH_5','Cond_3','Cond_5','pCO2_1','pCO2_2','pCO2_3','pCO2_4','pCO2_5','K_1','K_2','K_3','K_4','K_5','CO2_void','Corrected_Cond_3']
        res_all = pd.DataFrame(res,index = cols).T
        end = len(load_range)
        for i in range(1,end):
            self.q_area = I_dens / (load_range[i] * K_1 * ab.Fa)
            res = self.all_state(I_dens)
            resp = pd.DataFrame(res,index = cols).T
            frames = [res_all, resp]
            res_all = pd.concat(frames, ignore_index=True)
        return res_all
    
    def show(self,file_name):
        print(self.sol_1.name.center(80,'+'))
        print(*self.sol_1.state(),sep ='\n')
        print('characteristics'.center(80,'+'))
        print(*self.sol_1.characteristics(),sep= '\n')
        self.sol_1.save_excell(file_name)
        return

# t_start = time.time()
name = 'this should be ok '
cell_char = [2.5*10**-5,5*10**-3,7.5*10**-4,5*10**-3,1.075*10**-2,2.67,1]
anions = [ab.Compound('St',['H2SO4', 'HSO4', 'SO4'])]
anions_in = [0.05]
ed_test = ed_cell_ideal(1.1,1,anions,anions_in,cell_char)
ed_test.show(name)
print(*ed_test.all_state(0.75),sep='\n')
res_all = ed_test.trace()
res_all.to_csv(name+'.csv')
#print(res_all)
# t_middle = time.time()
#fp.plot_all([name,'Load_R [1]'],res_all,'Load_R',n_col=4,filesave =1)    
# t_finished = time.time()
# print('timing'.center(80,'+'))
# print(t_middle-t_start, sep ='\n')
# print(t_finished-t_middle, sep ='\n')