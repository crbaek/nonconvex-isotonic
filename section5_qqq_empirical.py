import numpy as np, pandas as pd, os, math
from scipy import sparse
from scipy.optimize import minimize
from scipy.stats import binomtest
from pandas.tseries.offsets import DateOffset

DATA='/mnt/data/tickhistoryintradaysummaries_5m_0x09d3fd799dbbe927_qqq-o_20200101000000000000000_20260426235959999999999.csv'
OUT='/mnt/data/section5_twosided_output'; os.makedirs(OUT,exist_ok=True)
TAUS=np.array([.95,.975,.99]); PLUS=np.array([.005,.0075,.010,.0125,.015]); MINUS=np.array([.01,.02,.03,.04,.05,.06]); CH=.80; WH=.06; CU=.90; WU=.18
DATE0=pd.Timestamp('2020-01-02'); DATE1=pd.Timestamp('2025-12-31'); FIRST=pd.Timestamp('2022-04-01')

# deterministic scalar response lookup; exact enough for empirical work, with local refinement
from functools import lru_cache
@lru_cache(None)
def response_table(kind,eta):
    u=np.linspace(0,1,10001); z=np.linspace(0,1,4001)
    # chunk to bound memory
    ans=np.empty_like(u)
    bump=-.5*z*z + (eta*np.exp(-((z-CH)**2)/(2*WH**2)) if kind=='plus' else -eta*np.exp(-((z-CU)**2)/(2*WU**2)))
    for a in range(0,len(u),250):
        uu=u[a:a+250,None]; vals=uu*z[None,:]+bump[None,:]
        # np.argmax gives first; maximal tie extremely rare; reverse for maximal tie
        idx=(vals.shape[1]-1)-np.argmax(vals[:,::-1],axis=1)
        ans[a:a+len(idx)]=z[idx]
    return u,ans

def T_resp(u,kind,eta):
    if kind=='zero' or eta==0: return np.clip(u,0,1)
    ug,zg=response_table(kind,float(eta)); return np.interp(np.clip(u,0,1),ug,zg)

def edges(G,K=3):
    def ix(k,m,q): return q*G*G+k*G+m
    E=[]
    for q in range(K):
      for k in range(G):
       for m in range(G):
        i=ix(k,m,q)
        if k<G-1:E.append((i,ix(k+1,m,q)))
        if m<G-1:E.append((i,ix(k,m+1,q)))
    for q in range(K-1):
      for k in range(G):
       for m in range(G): E.append((ix(k,m,q),ix(k,m,q+1)))
    return E

def iso_project(g,w,G):
    # Hildreth/Dykstra coordinate algorithm for weighted order-cone projection
    E=edges(G); x=g.astype(float).copy(); w=np.maximum(w.astype(float),1e-6); lam=np.zeros(len(E))
    for sweep in range(20000):
        maxv=0.0; maxchg=0.0
        for e,(i,j) in enumerate(E):
            v=x[i]-x[j]; maxv=max(maxv,v)
            den=1.0/w[i]+1.0/w[j]
            new=max(0.0, lam[e] + v/den)
            dl=new-lam[e]
            if dl!=0.0:
                x[i]-=dl/w[i]; x[j]+=dl/w[j]; lam[e]=new; maxchg=max(maxchg,abs(dl))
        if maxv<1e-9 and maxchg<1e-9: break
    viol=max(x[i]-x[j] for i,j in E)
    if viol>2e-6: raise RuntimeError(f'projection violation {viol} after {sweep+1} sweeps')
    return np.clip(x,0,1)

def q8(x,p): return np.quantile(x,p,method='median_unbiased') if len(x) else np.nan

def breaks(x,G):
    b=np.quantile(x,np.linspace(0,1,G+1),method='median_unbiased');
    # strict tiny jitter
    for i in range(1,len(b)):
      if b[i]<=b[i-1]: b[i]=np.nextafter(b[i-1],np.inf)
    b[0]=-np.inf;b[-1]=np.inf; return b

def bins(x,b): return np.clip(np.searchsorted(b,x,side='right')-1,0,len(b)-2)

def signals(d,G,rb,db,scale,fit=True):
    r=bins(d.rv1h.values,rb); m=bins(d.dd1h.values,db); y=np.clip(d.loss.values/scale,0,1)
    n=G*G*len(TAUS); g=np.empty(n); w=np.empty(n); stress=np.zeros(n,bool)
    glob=[q8(y,t) for t in TAUS]
    for q,t in enumerate(TAUS):
      for k in range(G):
       for j in range(G):
        ix=q*G*G+k*G+j; yy=y[(r==k)&(m==j)]; g[ix]=q8(yy,t) if len(yy) else glob[q]; w[ix]=len(yy) if len(yy) else (1e-6 if fit else 0); stress[ix]=(k>=math.floor(.75*G) and j>=math.floor(.75*G))
    return g,w,stress

def obs_loss(d,fit,G,rb,db,scale):
    r=bins(d.rv1h.values,rb); m=bins(d.dd1h.values,db); y=d.loss.values
    out=[]
    for q,t in enumerate(TAUS):
      pred=fit[q*G*G+r*G+m]*scale; e=y-pred; pin=np.mean(np.maximum(t*e,(t-1)*e)); hit=np.mean(y>pred)
      sm=(r>=math.floor(.75*G))&(m>=math.floor(.75*G)); esp=e[sm]; spin=np.mean(np.maximum(t*esp,(t-1)*esp)) if sm.any() else np.nan; shit=np.mean(y[sm]>pred[sm]) if sm.any() else np.nan
      out.append((t,pin,spin,hit,shit,int(len(y)),int(sm.sum())))
    return out

def prep(h):
    use=['Date-Time','Last','Close Bid','Close Ask','Close Mid Price']
    d=pd.read_csv(DATA,usecols=lambda c:c in use)
    ts=pd.to_datetime(d['Date-Time'],utc=True,errors='coerce').dt.tz_convert('America/New_York'); d['ts']=ts
    d=d[(ts.dt.tz_localize(None)>=DATE0)&(ts.dt.tz_localize(None)<DATE1+pd.Timedelta(days=1))]
    tm=d.ts.dt.strftime('%H:%M:%S'); d=d[(tm>='09:30:00')&(tm<='16:00:00')].copy()
    last=pd.to_numeric(d.get('Last'),errors='coerce'); bid=pd.to_numeric(d.get('Close Bid'),errors='coerce'); ask=pd.to_numeric(d.get('Close Ask'),errors='coerce'); mid=(bid+ask)/2
    if 'Close Mid Price' in d: cm=pd.to_numeric(d['Close Mid Price'],errors='coerce')
    else: cm=pd.Series(np.nan,index=d.index)
    d['price']=last.where(last>0,cm.where(cm>0,mid)); d=d[d.price>0].sort_values('ts').drop_duplicates('ts',keep='last')
    d['date']=d.ts.dt.date; d['logp']=np.log(d.price)
    # group computations, require exact spacing
    outs=[]; hb=h//5
    for day,x in d.groupby('date',sort=True):
      x=x.sort_values('ts').copy(); lp=x.logp.to_numpy(); tt=x.ts.to_numpy(); N=len(x)
      r=np.r_[np.nan,np.diff(lp)]; rv=np.full(N,np.nan); dd=np.full(N,np.nan); loss=np.full(N,np.nan)
      for i in range(12,N):
        # Match the submitted rolling-state construction: 12 most recent 5-minute returns
        # and 13 price observations within each trading day. Forward horizons below
        # are still required to span exactly h clock minutes.
        if np.all(np.isfinite(r[i-11:i+1])):
          rv[i]=np.sqrt(np.sum(r[i-11:i+1]**2))
        seg=lp[i-12:i+1]; dd[i]=max(0,max(seg[a]-np.min(seg[a:]) for a in range(len(seg))))
      if N>hb:
        for i in range(N-hb):
          if (x.ts.iloc[i+hb]-x.ts.iloc[i]).total_seconds()==h*60: loss[i]=-(lp[i+hb]-lp[i])
      x['rv1h']=rv;x['dd1h']=dd;x['loss']=loss; outs.append(x[['ts','date','rv1h','dd1h','loss']])
    z=pd.concat(outs,ignore_index=True).dropna(); return z

rows=[]; selrows=[]
for h in [int(os.environ.get('H','30'))]:
  d=prep(h); print('prepared',h,len(d),flush=True)
  starts=pd.date_range(FIRST,pd.Timestamp('2025-11-01'),freq='2MS')
  for fold,test0 in enumerate(starts,1):
    test1=test0+DateOffset(months=2)-pd.Timedelta(days=1); val0=test0-DateOffset(months=2); val1=test0-pd.Timedelta(days=1); tr0=val0-DateOffset(months=24); tr1=val0-pd.Timedelta(days=1)
    dates=pd.to_datetime(d['date']); tr=d[(dates>=tr0)&(dates<=tr1)]; va=d[(dates>=val0)&(dates<=val1)]; te=d[(dates>=test0)&(dates<=test1)]
    for G in [8,10]:
      rb=breaks(tr.rv1h.values,G); db=breaks(tr.dd1h.values,G); scale=max(q8(tr.loss.values,.999),1e-8)
      gt,wt,stress=signals(tr,G,rb,db,scale,True); proj=iso_project(gt,wt,G)
      gv,wv,_=signals(va,G,rb,db,scale,False)
      # validation weighted stress surface mse; test-independent
      candidates=[('zero',0.0)]+[('plus',float(e)) for e in PLUS]+[('minus',float(e)) for e in MINUS]
      vals=[]; mask=stress&(wv>0)
      for kind,eta0 in candidates:
        fit=T_resp(proj,kind,eta0); mse=np.sum(wv[mask]*(fit[mask]-gv[mask])**2)/np.sum(wv[mask]); vals.append(mse)
      vals=np.array(vals); conv=vals[0]
      # Hierarchical response selection: first identify direction using one fixed
      # representative strength per branch, then tune eta only within that branch.
      rep_idx={'plus':candidates.index(('plus',0.010)),'minus':candidates.index(('minus',0.030))}
      rep_scores={'zero':conv,'plus':vals[rep_idx['plus']],'minus':vals[rep_idx['minus']]}
      branch=min(rep_scores,key=rep_scores.get)
      if branch!='zero' and rep_scores[branch] > conv*(1-.005): branch='zero'
      if branch=='zero':
        si=0
      else:
        inds=np.array([i for i,c in enumerate(candidates) if c[0]==branch],dtype=int)
        best=vals[inds].min(); close=inds[vals[inds]<=best*(1+.005)]; si=int(close[0])
        if vals[si] > conv*(1-.005): si=0
      kind,eta=candidates[si]; selrows.append(dict(horizon=h,fold=fold,G=G,branch=kind,eta=eta,val_conv=conv,val_sel=vals[si],direction_branch=branch))
      # refit projection on train+validation, with grid/scale fixed from train (no test leakage)
      tv=pd.concat([tr,va]); gr,wr,_=signals(tv,G,rb,db,scale,True); pr=iso_project(gr,wr,G)
      fits={'convex':T_resp(pr,'zero',0),'two_sided_selected':T_resp(pr,kind,eta)}
      gtest,wtest,st=signals(te,G,rb,db,scale,False)
      for meth,fit in fits.items():
        mask=wtest>0; sm=st&mask
        gm=np.sum(wtest[mask]*(fit[mask]-gtest[mask])**2)/np.sum(wtest[mask]); smse=np.sum(wtest[sm]*(fit[sm]-gtest[sm])**2)/np.sum(wtest[sm])
        for tau,pin,spin,hit,shit,nobs,nstress in obs_loss(te,fit,G,rb,db,scale):
          rows.append(dict(horizon=h,fold=fold,G=G,method=meth,branch=kind,eta=eta,tau=tau,grid_mse=gm,stress_mse=smse,pinball=pin,stress_pinball=spin,hit=hit,stress_hit=shit,nobs=nobs,nstress=nstress))
    print('h',h,'fold',fold,'done',flush=True)

raw=pd.DataFrame(rows); sel=pd.DataFrame(selrows); raw.to_csv(OUT+f'/raw_h{h}.csv',index=False);sel.to_csv(OUT+f'/selection_h{h}.csv',index=False)
# paired gains per fold metric, q-specific for pinball
pairs=[]
for keys,x in raw.groupby(['horizon','fold','G','tau']):
 c=x[x.method=='convex'].iloc[0]; s=x[x.method=='two_sided_selected'].iloc[0]
 rec=dict(zip(['horizon','fold','G','tau'],keys)); rec['eta']=s.eta; rec['branch']=s.branch
 for met in ['grid_mse','stress_mse','pinball','stress_pinball']:
  rec[met+'_gain']=100*(c[met]-s[met])/c[met] if c[met]!=0 else np.nan
 rec['hit_conv']=c.hit;rec['hit_new']=s.hit;rec['stress_hit_conv']=c.stress_hit;rec['stress_hit_new']=s.stress_hit
 pairs.append(rec)
pair=pd.DataFrame(pairs); pair.to_csv(OUT+f'/paired_h{h}.csv',index=False)
# aggregate
agg=[]
for keys,x in pair.groupby(['horizon','G','tau']):
 rec=dict(zip(['horizon','G','tau'],keys)); rec['nfold']=len(x); rec['p_nonconvex']=np.mean(x.branch!='zero'); rec['p_plus']=np.mean(x.branch=='plus'); rec['p_minus']=np.mean(x.branch=='minus')
 for met in ['grid_mse_gain','stress_mse_gain','pinball_gain','stress_pinball_gain']:
  rec[met]=x[met].mean();rec[met+'_se']=x[met].std(ddof=1)/np.sqrt(len(x))
 rec['hit_conv']=x.hit_conv.mean();rec['hit_new']=x.hit_new.mean();rec['stress_hit_conv']=x.stress_hit_conv.mean();rec['stress_hit_new']=x.stress_hit_new.mean(); agg.append(rec)
pd.DataFrame(agg).to_csv(OUT+f'/summary_h{h}.csv',index=False)
print(pd.DataFrame(agg).to_string(index=False),flush=True)
