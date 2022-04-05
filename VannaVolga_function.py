import BlackScholes
from scipy.stats import norm
from scipy import optimize
from math import pow,log,exp,sqrt
from pandas import read_csv


def VolK(v_atm,BF,RR,CallPut):
    return v_atm +BF +BlackScholes.cp(CallPut.upper())*RR/2

def GetStrikeFromDelta(S,fwd,v,expiry_date,value_date,CallPut,delta):
    c_p = BlackScholes.cp(CallPut.upper())
    f = (S+fwd)
    t = BlackScholes.T(expiry_date,value_date)

    if v > 1.5:
        v = v/100
    else:
        v = v
    return f * exp(-c_p*norm.ppf(delta)*v*sqrt(t)+((pow(v,2))/2)*t)


def GetStrikeFromDeltaPA(S,fwd,v,expiry_date,value_Date,CallPut,delta):

    c_p = BlackScholes.cp(CallPut.upper())
    f= (S+fwd)
    t = BlackScholes.T(expiry_date,value_Date)

    K_unadjust = GetStrikeFromDelta(S,fwd,v,expiry_date,value_Date,CallPut,delta)
    get_k= lambda K : abs(BlackScholes.BlackScholes(S,K,fwd,v,expiry_date,value_Date,CallPut,'BUY',1,True).delta()) - 0.25
    return optimize.newton(get_k,K_unadjust)

def Get_DN_Strike(S,fwd,v,expiry_date,value_date,deltaBase=True):

    guess_k = (S+fwd)*exp(-norm.ppf(0.5)*(v/100)*sqrt(BlackScholes.T(expiry_date,value_date))+(pow((v/100),2))/2)
    def opt_call(S,K,fwd,v,expiry_date,value_date,deltaBase):
        return BlackScholes.BlackScholes(S,K,fwd,v,expiry_date,value_date,'CALL','BUY',1,deltaBase).delta()
    def opt_put(S,K,fwd,v,expiry_date,value_date,deltaBase):
        return BlackScholes.BlackScholes(S,K,fwd,v,expiry_date,value_date,'PUT','BUY',1,deltaBase).delta()
    get_k = lambda K : opt_call(S,K,fwd,v,expiry_date,value_date,deltaBase) + opt_put(S,K,fwd,v,expiry_date,value_date,deltaBase)
    return optimize.newton(get_k,guess_k)

def convert_bfly(v,BF,RR):
    return BF+(-0.02*RR)+(0.0005778*abs(RR/BF))+(-0.17*RR/v)

def load_data(file_location,file_name):
    location = file_location
    fname = file_name
    df = read_csv(fr'{location}/{fname}')
    df['Expiry'] = df['expiry']
    df= df.set_index('expiry')
    return df

def interploate_data(expiry_date,value_date,df):
    days_to_expiry = BlackScholes.T(expiry_date,value_date)

    df['Time_to_Expiry'] = df['Expiry'].apply(lambda x:BlackScholes.T(x,value_date))
    expiries = df['Time_to_Expiry'].to_numpy()
    market_data = df.to_numpy()

    idx = 0
    for date in expiries:
        if days_to_expiry > date:
            if days_to_expiry < expiries[idx+1]:
                w1 = (expiries[idx+1] - days_to_expiry)/(expiries[idx+1]-expiries[idx])
                w2 = 1-w1
                v1 = market_data[idx,1]**2
                v2 = market_data[idx+1,1]**2
                v = sqrt(v2*w1 + v1*w2)
                fwd = (market_data[(idx+1),0]*w1 + market_data[idx,0]*w2)
                rr_25d = (market_data[(idx+1),2]*w1+market_data[idx,2]*w2)
                rr_10d = (market_data[(idx + 1),3] * w1 + market_data[idx,3] * w2)
                bf_25d = (market_data[(idx + 1),4] * w1 + market_data[idx,4] * w2)
                bf_10d = (market_data[(idx + 1),5] * w1 + market_data[idx,5] * w2)

                break
            elif days_to_expiry > expiries[-1]:
                print('Date is out of range')
                break
            else:
                idx+=1
    return v,fwd,rr_25d,rr_10d,bf_25d,bf_10d
