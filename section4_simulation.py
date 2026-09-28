import numpy as np, pandas as pd, time, os
from scipy.optimize import minimize, LinearConstraint, minimize_scalar, brentq
from scipy.stats import norm, t as tdist
from numba import njit

B=int(os.getenv('SIM_B','100')); BASE=20260918
K=L=12; alphas=np.array([.95,.99]); ntrain,nval,ntest=60,120,200
cH,wH=.80,.06; etas=np.array([0,.005,.0075,.010,.0125,.015])
scenarios=[('smooth',0),('weak_transition',.010),('pronounced_transition',.015)]
errors=['gaussian','t5','t3']
if os.getenv('ONLY_SIGNAL'): scenarios=[x for x in scenarios if x[0]==os.getenv('ONLY_SIGNAL')]
if os.getenv('ONLY_ERROR'): errors=[os.getenv('ONLY_ERROR')]
idx=np.arange(K*L*2).reshape((K,L,2),order='F')
edges=[]
for tt in range(2):
 for k in range(K-1):
  for m in range(L): edges.append((idx[k,m,tt],idx[k+1,m,tt]))
 for k in range(K):
  for m in range(L-1): edges.append((idx[k,m,tt],idx[k,m+1,tt]))
for k in range(K):
 for m in range(L): edges.append((idx[k,m,0],idx[k,m,1]))
edges=np.array(edges,dtype=int); n=K*L*2
A=np.zeros((len(edges),n)); rr=np.arange(len(edges)); A[rr,edges[:,1]]=1; A[rr,edges[:,0]]=-1
LC=LinearConstraint(A,0,np.inf)

# Fast independent isotonic projection via Dykstra's algorithm on the order halfspaces.
# This solves the same Euclidean projection without using Gurobi's quadratic-objective convention.
@njit(cache=True)
def _dykstra(g, edges, tol=1e-9, max_cycles=5000):
    x=g.copy(); p=np.zeros(edges.shape[0])
    for cyc in range(max_cycles):
        md=0.0
        for r in range(edges.shape[0]):
            i=edges[r,0]; j=edges[r,1]; oldi=x[i]; oldj=x[j]; a=p[r]
            yi=oldi+a; yj=oldj-a; v=yi-yj
            if v>0.0:
                d=0.5*v; x[i]=yi-d; x[j]=yj+d; p[r]=d
            else:
                x[i]=yi; x[j]=yj; p[r]=0.0
            di=abs(x[i]-oldi); dj=abs(x[j]-oldj)
            if di>md: md=di
            if dj>md: md=dj
        if md<tol: return x,cyc+1
    return x,max_cycles

def iso(g):
    y,cyc=_dykstra(np.asarray(g,dtype=np.float64),edges)
    viol=np.max(y[edges[:,0]]-y[edges[:,1]])
    if viol>2e-7: raise RuntimeError(f'order violation {viol}')
    return y

def scalar_obj(z,u,e):
    return u*z-.5*z*z+e*np.exp(-((z-cH)**2)/(2*wH*wH))

def T_one(u,e):
    if e==0: return float(np.clip(u,0,1))
    # Find all stationary points from derivative sign changes on a fine grid,
    # then compare every stationary point with both boundaries: global solve.
    z=np.linspace(0,1,801)
    ex=np.exp(-((z-cH)**2)/(2*wH*wH))
    d=u-z-e*(z-cH)/(wH*wH)*ex
    roots=[]
    exact=np.where(np.abs(d)<1e-12)[0]
    roots.extend(z[exact].tolist())
    ch=np.where(d[:-1]*d[1:]<0)[0]
    def der(x): return u-x-e*(x-cH)/(wH*wH)*np.exp(-((x-cH)**2)/(2*wH*wH))
    for j in ch:
        roots.append(brentq(der,z[j],z[j+1],xtol=1e-13))
    cand=np.array([0.,1.]+roots)
    vals=np.array([scalar_obj(x,u,e) for x in cand])
    mx=vals.max()
    # maximal maximizer tie rule
    return float(cand[vals>=mx-1e-10].max())

def T_vec(u,e):
    u=np.asarray(u); out=np.empty_like(u,dtype=float)
    # pooled values repeat heavily; solve each unique value once
    ur=np.round(u,12); vals=np.unique(ur)
    mp={v:T_one(float(v),e) for v in vals}
    for j,v in enumerate(ur): out[j]=mp[v]
    return out

def base_u():
    a=np.zeros((K,L,2))
    for tt in range(2):
      for k in range(K):
       for m in range(L): a[k,m,tt]=.16+.30*k/(K-1)+.20*m/(L-1)+[0,.10][tt]
    return a

def target(e):
    u=base_u(); return T_vec(u.ravel(order='F'),e).reshape(u.shape,order='F')

def qerr(p,error):
    if error=='gaussian': return norm.ppf(p)
    return tdist.ppf(p,df=int(error[1:]))

def calibrate(tar,error):
    q1,q2=qerr(alphas,error)
    sig=(tar[:,:,1]-tar[:,:,0])/(q2-q1)
    mu=tar[:,:,0]-sig*q1
    if np.any(sig<=0): raise RuntimeError('nonpositive scale')
    return mu,sig

def losses(mu,sig,N,error,rng):
    if error=='gaussian': e=rng.normal(size=(K,L,N))
    else: e=rng.standard_t(int(error[1:]),size=(K,L,N))
    return mu[:,:,None]+sig[:,:,None]*e

def sampleq(x):
    q=np.quantile(x,alphas,axis=2,method='linear') # shape 2,K,L
    out=np.moveaxis(q,0,2)
    out[:,:,1]=np.maximum(out[:,:,1],out[:,:,0])
    return out

def stressmask(tar): return tar[:,:,1]>=cH-wH

def stress_mse(a,b,mask):
    aa=np.concatenate([a[:,:,0][mask],a[:,:,1][mask]])
    bb=np.concatenate([b[:,:,0][mask],b[:,:,1][mask]])
    return np.mean((aa-bb)**2)
def pinball(x,fit):
    out=[]
    for tt,a in enumerate(alphas):
      u=x-fit[:,:,tt,None]; out.append(np.mean(u*(a-(u<0))))
    return out

def runone(rep,signal,error,edgp):
    sid=(['smooth','weak_transition','pronounced_transition'].index(signal)+1)*10+(['gaussian','t5','t3'].index(error)+1)
    # R set.seed cannot be replicated exactly; independent random stream is intentional.
    rng=np.random.default_rng(BASE+10000*sid+rep)
    tar=target(edgp); mu,sig=calibrate(tar,error)
    tr=losses(mu,sig,ntrain,error,rng); va=losses(mu,sig,nval,error,rng); te=losses(mu,sig,ntest,error,rng)
    qtr=sampleq(tr); qva=sampleq(va)
    g=np.clip(qtr.ravel(order='F'),0,1); gs=iso(g); mask=stressmask(tar)
    cand=[]
    fits={}
    for e in etas:
      xv=T_vec(gs,e).reshape((K,L,2),order='F'); fits[e]=xv
      cand.append((e,stress_mse(xv,qva,mask)))
    base=cand[0][1]; best=min(x[1] for x in cand)
    near=[x for x in cand if x[1]<=best*1.005]; es=min(near,key=lambda x:x[0])[0]
    sval=dict(cand)[es]
    if es>0 and sval>base*(1-.005): es=0.
    rows=[]
    for meth,e in [('convex',0.),('stress_selected',float(es))]:
      x=fits[e]; ql=pinball(te,x)
      rows.append(dict(rep=rep,signal=signal,error=error,eta_dgp=edgp,method=meth,eta=e,
        grid_mse=np.mean((x-tar)**2),stress_mse=stress_mse(x,tar,mask),qloss_095=ql[0],qloss_099=ql[1],
        max_order_violation=max(0,np.max(x.ravel(order='F')[edges[:,0]]-x.ravel(order='F')[edges[:,1]]))))
    return rows

# audits
ordered=np.linspace(.1,.9,n)
y=iso(ordered)
assert np.max(np.abs(y-ordered))<1e-7
# factor-of-two diagnostic: unconstrained/feasible identity is preserved by construction
print('Preflight identity max error',np.max(np.abs(y-ordered)))

rows=[]; t0=time.time()
for signal,edgp in scenarios:
 for error in errors:
  for rep in range(1,B+1):
   rows.extend(runone(rep,signal,error,edgp))
  print(signal,error,'done',round(time.time()-t0,1),'sec',flush=True)
df=pd.DataFrame(rows); prefix=os.getenv('OUT_PREFIX','section4_python'); df.to_csv('/mnt/data/'+prefix+'_raw_metrics.csv',index=False)
# paired
w=df.pivot(index=['rep','signal','error','eta_dgp'],columns='method',values=['grid_mse','stress_mse','qloss_095','qloss_099','eta'])
w.columns=['_'.join(c) for c in w.columns]; w=w.reset_index()
for v in ['grid_mse','stress_mse','qloss_095','qloss_099']:
 w[v+'_gain_pct']=100*(w[v+'_convex']-w[v+'_stress_selected'])/w[v+'_convex']
w.to_csv('/mnt/data/'+prefix+'_paired_gains.csv',index=False)
summary=w.groupby(['signal','error','eta_dgp']).agg(
 grid_mse_gain_pct=('grid_mse_gain_pct','mean'),grid_mse_gain_se=('grid_mse_gain_pct',lambda x:x.std(ddof=1)/np.sqrt(len(x))),
 stress_mse_gain_pct=('stress_mse_gain_pct','mean'),stress_mse_gain_se=('stress_mse_gain_pct',lambda x:x.std(ddof=1)/np.sqrt(len(x))),
 qloss_095_gain_pct=('qloss_095_gain_pct','mean'),qloss_099_gain_pct=('qloss_099_gain_pct','mean'),
 positive_eta_rate=('eta_stress_selected',lambda x:np.mean(x>0)),mean_selected_eta=('eta_stress_selected','mean')).reset_index()
summary.to_csv('/mnt/data/'+prefix+'_gain_summary.csv',index=False)
print(summary.to_string(index=False))
print('TOTAL SEC',time.time()-t0)
