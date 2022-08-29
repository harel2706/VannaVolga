import numpy as np
from math import log,sqrt,exp,pi
import datetime,time
from scipy.stats import norm
from datetime import timedelta

def T(expiry_date,value_date):
    date2 = datetime.datetime.strptime(expiry_date,'%Y-%m-%d')
    date1= datetime.datetime.strptime(value_date,'%Y-%m-%d')
    return ((date2-date1)/timedelta(days=365))

def cp(CallPut):
    return {'CALL':1,'PUT':-1}[CallPut]


def pos(BuySell):
    return {'BUY':1,'SELL':-1}[BuySell]


class BlackScholes:

    def __init__(self,S,K,fwd,v,expiry_date,value_date,CallPut,BuySell,Notional,deltaBase=True):
        self.S = S
        self.K = K
        self.f = (S+fwd)
        self.v = v/100
        self.expiry_date = expiry_date
        self.value_date = value_date
        self.t = T(expiry_date,value_date)
        self.CallPut = cp(CallPut.upper())
        self.BuySell = pos(BuySell.upper())
        self.Notional = Notional
        self.deltaBase = deltaBase

    def d1(self):
        f,k,v,t = self.f,self.K,self.v,self.t
        return (log(f/k)+(pow(v,2)*t)/2)/(sqrt(t)*v)

    def d2(self):
        v,t = self.v,self.t
        return BlackScholes.d1(self) - v*sqrt(t)

    def price(self):
        f,k,v,t,c_p,bs,n = self.f,self.K,self.v,self.t,self.CallPut,self.BuySell,self.Notional
        return bs*c_p*(f*norm.cdf(c_p*BlackScholes.d1(self))-k*norm.cdf(c_p*BlackScholes.d2(self)))*n/self.S

    def delta(self):
        bs , n ,c_p =self.BuySell,self.Notional ,self.CallPut

        delta = np.where(self.deltaBase==True,
                         (bs *(c_p* norm.cdf(c_p* BlackScholes.d1(self)) - BlackScholes.price(self))*n),
                         (bs *c_p * norm.cdf(c_p* BlackScholes.d1(self))*n))

        return delta

    def vega(self):
        s,f,t,bs,n = self.S,self.f,self.t,self.BuySell,self.Notional
        return bs * ((norm.pdf(BlackScholes.d1(self)))*f*sqrt(t)*0.01/s)*n

    def theta(self):
        s,f,t,bs,n = self.S,self.f,self.t,self.BuySell,self.Notional
        return round(bs * (((-1/365*f*norm.pdf(BlackScholes.d1(self)))/2*sqrt(t))/s)*n,2)

    #2nd order greeks

    def gamma(self):
        f,v,t,bs,n = self.f, self.v,self.t, self.BuySell, self.Notional
        return round(bs * ((norm.pdf(BlackScholes.d1(self))*f*0.01)/(f*v*sqrt(t)))*n,2)

    def volga(self):
        v = self.v
        return ((self.BuySell**2)*BlackScholes.vega(self)*BlackScholes.d1(self)*BlackScholes.d2(self))

    def vanna(self):
        s ,v ,t ,bs  = self.S,self.v,self.t,self.BuySell
        return (pow(bs,2))*(BlackScholes.vega(self)/s)*(1-(BlackScholes.d1(self)/(v * t)))

    def dVdK(self):
        f, k, v ,t  = self.f , self.K , self.v ,self.t

        sqr_2pi = sqrt(2*pi*(pow(v,2))*t)
        d1 = pow((log(k/f)+((pow(v,2))/2)*t),2)
        return (1/k)*(1/sqr_2pi)*exp(-(1/(2*(v**2)*t))*(d1))

    def vega_volga_risk(self):
        vega = BlackScholes.vega(self)
        volga = BlackScholes.volga(self)

        return (vega+volga)

