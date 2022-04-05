import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns

import VannaVolga_ImpliedVol
from BlackScholes import BlackScholes , T
from VannaVolga_function import GetStrikeFromDeltaPA , GetStrikeFromDelta , VolK , load_data,interploate_data
from VannaVolga_ImpliedVol import VannaVolga
sns.set()

import time


t0 = time.time()

delta_call_put =[]
delta_range = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49, 0.49, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2,
                   0.15, 0.1, 0.05]
delta_names = ['5d Put', '10d Put', '15d Put', '20d Put', '25d Put', '30d Put', '35d Put', '40d Put', '45d Put',
                   '50d Put',
                   '50d Call', '45d Call', '40d Call', '35d Call', '30d Call', '25d Call', '20d Call', '15d Call',
                   '10d Call', '5d Call']
call_put = ['PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'PUT', 'CALL', 'CALL', 'CALL',
                'CALL', 'CALL', 'CALL', 'CALL', 'CALL', 'CALL', 'CALL']


delta_strike = [[delta,cp] for delta,cp in zip(delta_range,call_put)]
delta_call_put.append(delta_strike)
delta_call_put = delta_call_put[0]


def plot_smile(S,fwd,v,RR,BF,expiry_date,value_date,delta,deltaBase=True,atm_conv='DN'):
    vv_vol = []

    for strike in delta_call_put:
        K = GetStrikeFromDelta(S, fwd, v, expiry_date, value_date, strike[1], strike[0])
        vol = VannaVolga(S,K,fwd,expiry_date,value_date,v,RR,BF,'CALL',0.25,deltaBase,atm_conv).GetImpliedVol()
        vv_vol.append(vol)

    plt.plot(vv_vol)
    plt.xticks(np.arange(0,len(delta_names)),delta_names,rotation=45,fontsize=8)
    plt.show()

plot_smile(114.29,-0.170,6.26,-0.77,0.225,'2022-03-22','2021-12-22',0.25,True,'DN')

def Get_Vol_from_file(file_location,file_name,ref_spot,strike,expiry,value_date,fwd_div=10000,deltaBase=True,atm_conv='DN'):
    df = load_data(file_location,file_name)
    if expiry not in df.index:
        interp_data = interploate_data(expiry, value_date, df)
        v = interp_data[0]
        fwd = interp_data[1]/fwd_div
    else:
        market_data = df.loc[expiry,:]
        fwd = market_data['fwd']/fwd_div
        v = market_data['atm_vol']

    bs_delta = BlackScholes(S=ref_spot,K=strike,v=v,fwd=fwd,
                            expiry_date=expiry_date,value_date=value_date,
                            CallPut='CALL',BuySell='BUY',Notional=1,deltaBase=deltaBase).delta()
    bs_delta= np.where(bs_delta>0.5,bs_delta-0.5,bs_delta)

    if expiry not in df.index:
        rr = np.where(bs_delta > 0.15,interp_data[2],interp_data[3])
        bf= np.where(bs_delta > 0.15,interp_data[4],interp_data[5])
    else:
        rr= np.where(bs_delta > 0.15 ,market_data['rr_25d'] ,market_data['rr_10d']  )
        bf = np.where(bs_delta > 0.15, market_data['bf_25d'],market_data['bf_10d'])
    delta = np.where(bs_delta > 0.15 , 0.25 , 0.10)
    vol = VannaVolga_ImpliedVol.VannaVolga(S=ref_spot,K=strike,fwd=fwd,expiry_date=expiry,value_date=value_date,
                                           v_atm=v,RR=rr,BF=bf,CallPut='CALL',delta=delta,deltaBase=deltaBase,atm_conv=atm_conv).GetImpliedVol()

    print(vol)

def plot_surface(S,fwd_curve,v_atm,RR_term,BF_term,expiries,value_date,delta,deltaBase=True,atm_conv='DN'):

    market_data=[]
    vol_surface =[]

    data = [[expiry,fwd,vol,rr,bf] for expiry,fwd,vol,rr,bf in zip(expiries,fwd_curve,v_atm,RR_term,BF_term)]
    market_data.append(data)
    market_data = market_data[0]

    vol = [[VannaVolga(S,
                          GetStrikeFromDelta(S=S,fwd=term[1],v=term[2],
                                             expiry_date=term[0],value_date=value_date,CallPut=strike[1],
                                             delta=strike[0]),
                          fwd=term[1],expiry_date=term[0],value_date=value_date,
                          v_atm=term[2],RR=term[3],BF=term[4],CallPut='CALL',delta=delta,
                          deltaBase=deltaBase,atm_conv=atm_conv).GetImpliedVol() for strike in delta_call_put] for term in market_data]
    vol_surface.append(vol)


    vol_surface = np.array(vol_surface).reshape((len(expiries),len(delta_call_put)),order='f')
    df = pd.DataFrame(vol_surface,index=expiries,columns=delta_names)

    fig = go.Figure(data=[go.Surface(z=df.values,
                                     x=delta_names,
                                     y=expiries, colorscale='Viridis')])

    fig.update_layout(autosize=True,width=1000,height=1000,scene=dict(xaxis_title='delta strikes',
                                                      yaxis_title='expiries',
                                                      zaxis_title='volatility',))

    fig.show()

expiries = ['2021-12-29','2022-01-04','2022-01-22','2022-02-22','2022-03-22','2022-06-22','2022-09-22','2022-12-22']
fwd_curve = [0 ,-0.034,-0.038,-0.048,-0.07,-0.106,-.16,-0.210 ]
atm_curve = [4.8,5.07,5.67,6.02,6.25,6.48,6.59,6.71]
RR_curve = [-0.433,-0.66,-0.71,-0.76,-0.85,-0.888,-0.92,-0.92]
BF_curve = [0.16,0.18,0.2,0.203,0.23,0.255,0.285,0.305]

plot_surface(114.2,fwd_curve,atm_curve,RR_curve,BF_curve,expiries,'2021-12-20',0.25)

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


atm_vol_ref = 7.92
ref_spot = 3.1485
expiry_date='2022-04-21'
value_date  = '2022-01-21'

def smile_change(spot,ref_spot,fwd,ref_fwd,v_atm,v_atm_ref,rr,rr_ref,bf,bf_ref,expiry_date,value_date):

    fixed_strikes = [GetStrikeFromDelta(ref_spot, ref_fwd, atm_vol_ref, expiry_date, value_date, strike[1], strike[0]) for
                     strike in delta_call_put]

    ref_smile = [VannaVolga(ref_spot, strike, ref_fwd, expiry_date, value_date, v_atm_ref, rr_ref, bf_ref, 'CALL', 0.25).GetImpliedVol()
                 for strike in fixed_strikes]

    cur_smile = [VannaVolga(spot, strike, fwd, expiry_date, value_date, v_atm, rr, bf, 'CALL', 0.25).GetImpliedVol()
                 for strike in fixed_strikes]

    plt.plot(ref_smile)
    plt.plot(cur_smile, color='r', ls='--')

    plt.legend(['Ref. smile', 'Cur. smile'])
    plt.title('vol smile')
    plt.xticks(np.arange(0, len(fixed_strikes)), labels=np.around(fixed_strikes, 3), fontsize=8, rotation=45)
    plt.xlabel('fixed strikes')
    plt.ylabel('volatility')
    plt.show()


smile_change(1.04,1.0319,0,0,7.65,7.57,-0.985,-0.99,0.218,0.21,'2022-04-21','2022-01-21')

Get_Vol_from_file('C:/Users/user/Downloads','market_data.csv',1.1039,1.1270,'2022-06-04','2022-03-21',
                     fwd_div=10000,deltaBase=True,atm_conv='DN')

t1 = time.time()

print(f'total time elapsed :{t1- t0:.2f}')

