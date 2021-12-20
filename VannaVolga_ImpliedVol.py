import BlackScholes as bs
from BlackScholes import BlackScholes
import VannaVolga_function as vv
import numpy as np
import math
from scipy.stats import norm



class VannaVolga:

    def __init__(self,S,K,fwd,expiry_date,value_date,v_atm , RR , BF , CallPut,delta,deltaBase=True,atm_conv='DN'):
        self.S = S
        self.fwd = fwd
        self.f = (S+fwd)
        self.K = K
        self.expiry_date = expiry_date
        self.value_date = value_date
        self.t = bs.T(expiry_date,value_date)
        self.v_atm = v_atm
        self.v = v_atm/100
        self.RR = RR
        self.BF = BF
        self.CallPut = CallPut.upper()
        self.delta = delta
        self.deltaBase = deltaBase
        self.atm_conv = atm_conv

    # Strike Convetion

    def GetImpliedVol(self):
        s, fwd , f , K , t, expiry_date,value_date, v , RR , BF ,delta = \
            self.S, self.fwd ,self.f , self.K ,self. t ,self.expiry_date,self.value_date ,self.v_atm,self.RR , self.BF , self.delta

        def d1(s,fwd,k,v,t):
            return (np.log((s+fwd)/k)+((v**2)/2)*t)/(np.sqrt(t)*v)

        def d2(s,fwd,k,v,t):
            return d1(s,fwd,k,v,t) - np.sqrt(t)*v

        S1 = vv.VolK(v,BF,RR,'PUT')/100
        S2 = v/100
        S3 = vv.VolK(v,BF,RR,'CALL')/100

        if self.deltaBase==True:
            K1 = vv.GetStrikeFromDeltaPA(s,fwd,vv.VolK(v,BF,RR,'PUT'),expiry_date,value_date,'PUT',delta)
            if self.atm_conv=='DN':
                K2 = vv.Get_DN_Strike(s,fwd,v,expiry_date,value_date,deltaBase=True)
            elif self.atm_conv=='ATMF':
                K2 = f
            else:
                K2 = s
            K3 = vv.GetStrikeFromDeltaPA(s,fwd,vv.VolK(v,BF,RR,'CALL'),expiry_date,value_date,'CALL',delta)
        else:
            K1 = vv.GetStrikeFromDelta(s, fwd, vv.VolK(v, BF, RR, 'PUT'), expiry_date, value_date, 'PUT', delta)
            if self.atm_conv=='DN':
                K2 = vv.Get_DN_Strike(s,fwd,v,expiry_date,value_date,deltaBase=False)
            elif self.atm_conv=='ATMF':
                K2 = f
            else:
                K2 = f * np.exp((S2 ** 2) / 2 * t)
            K3 = vv.GetStrikeFromDelta(s, fwd, vv.VolK(v, BF, RR, 'CALL'), expiry_date, value_date, 'CALL', delta)

        y1 = (math.log(K2/K)*math.log(K3/K))/(math.log(K2/K1)*math.log(K3/K1))
        y2 = (math.log(K/K1)*math.log(K3/K))/(math.log(K2/K1)*math.log(K3/K2))
        y3 = (math.log(K/K1)*math.log(K/K2))/(math.log(K3/K1)*math.log(K3/K2))



        p = y1*S1 + y2 * S2 + y3* S3 - S2
        q = y1 * d1(s,fwd,K1,S2,t)*d2(s,fwd,K3,S2,t)*((S1-S2)**2) +\
            y3 * d1(s,fwd,K3,S2,t)*d2(s,fwd,K3,S2,t)*((S3-S2)**2)
        d1d2 = d1(s,fwd,K,S2,t)*d2(s,fwd,K,S2,t)
        return (S2+(-S2+np.sqrt(S2**2+d1d2*(2*S2*p+q)))/d1d2)*100

    def plot_surface(self):
        delta_call_put = []
        delta_range = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49, 0.49, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2,
                       0.15, 0.1, 0.05]
        delta_names = ['5d Put', '10d Put', '15d Put', '20d Put', '25d Put', '30d Put', '35d Put', '40d Put', '45d Put',
                       '50d Put',
                       '50d Call', '45d Call', '40d Call', '35d Call', '30d Call', '25d Call', '20d Call', '15d Call',
                       '10d Call', '5d Call']
        call_put = ['PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'CALL', 'CALL', 'CALL',
                    'CALL', 'CALL', 'CALL', 'CALL', 'CALL', 'CALL', 'CALL']



print(VannaVolga(113.53,112,-0.0511,'2022-01-20','2021-12-20',6.013,-0.853,0.210,'CALL',0.25,deltaBase=True,atm_conv='DN').GetImpliedVol())
