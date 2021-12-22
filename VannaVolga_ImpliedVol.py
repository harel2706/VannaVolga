import matplotlib.pyplot as plt
import seaborn as sns
import BlackScholes as bs
from BlackScholes import BlackScholes
import VannaVolga_function as vv
import numpy as np
import math
sns.set()


class VannaVolga:

    def __init__(self,S,K,fwd,expiry_date,value_date,v_atm , RR , BF , CallPut,delta,deltaBase=True,atm_conv='DN',strike_type='rate'):
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
        self.strike_type = strike_type

    # Strike Convetion

    def GetImpliedVol(self):
        s, fwd , f , K , t, expiry_date,value_date, v , RR , BF ,delta ,CallPut = \
            self.S, self.fwd ,self.f , self.K ,self. t ,self.expiry_date,self.value_date ,\
            self.v_atm,self.RR , self.BF , self.delta,self.CallPut

        def d1(s,fwd,k,v,t):
            return (np.log((s+fwd)/k)+((v**2)/2)*t)/(np.sqrt(t)*v)

        def d2(s,fwd,k,v,t):
            return d1(s,fwd,k,v,t) - np.sqrt(t)*v

        BF = vv.convert_bfly(v,BF,RR)

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

        if self.strike_type=='delta':
            if self.deltaBase==True:
                if CallPut=='CALL':
                    K = vv.GetStrikeFromDeltaPA(s,fwd,v,expiry_date,value_date,'CALL',K)
                else:
                    K = vv.GetStrikeFromDeltaPA(s,fwd,v,expiry_date,value_date,'PUT',K)
            else:
                if CallPut=='CALL':
                    K = vv.GetStrikeFromDelta(s,fwd,v,expiry_date,value_date,'CALL',K)
                else:
                    K = vv.GetStrikeFromDelta(s,fwd,v,expiry_date,value_date,'PUT',K)

        y1 = (math.log(K2/K)*math.log(K3/K))/(math.log(K2/K1)*math.log(K3/K1))
        y2 = (math.log(K/K1)*math.log(K3/K))/(math.log(K2/K1)*math.log(K3/K2))
        y3 = (math.log(K/K1)*math.log(K/K2))/(math.log(K3/K1)*math.log(K3/K2))



        p = y1*S1 + y2 * S2 + y3* S3 - S2
        q = y1 * d1(s,fwd,K1,S2,t)*d2(s,fwd,K3,S2,t)*((S1-S2)**2) +\
            y3 * d1(s,fwd,K3,S2,t)*d2(s,fwd,K3,S2,t)*((S3-S2)**2)
        d1d2 = d1(s,fwd,K,S2,t)*d2(s,fwd,K,S2,t)
        return (S2+(-S2+np.sqrt(S2**2+d1d2*(2*S2*p+q)))/d1d2)*100





print(VannaVolga(114.205,112,-0.1035,'2022-03-22','2021-12-22',6.26,-0.77,0.238,'CALL',0.25,deltaBase=True,atm_conv='DN',strike_type='rate').GetImpliedVol())


prob_surface,prob_flat_vol =[] , []
ref_spot = 12.475
fwd = 0.42
expiry_date = '2022-03-20'
value_date = '2021-12-20'
atm_vol= 74.67
RR = 10.28
BF = 3.27
strike_range = np.linspace(0.1,2.5)*ref_spot

for strike in strike_range:
    vol = VannaVolga(ref_spot,strike,fwd,expiry_date,value_date,atm_vol,RR,BF,'CALL',0.25,atm_conv='DN').GetImpliedVol()
    prob = BlackScholes(ref_spot,strike,fwd,vol,expiry_date,value_date,'CALL','BUY',1).dVdK()
    prob_flat = BlackScholes(ref_spot,strike,fwd,atm_vol,expiry_date,value_date,'CALL','BUY',1).dVdK()
    prob_surface.append(prob)
    prob_flat_vol.append(prob_flat)

plt.plot(prob_surface)
plt.plot(prob_flat_vol,color='r',ls='--')
plt.show()