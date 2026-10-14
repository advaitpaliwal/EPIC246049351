"""Copied from ../../scripts/fit.py (not modified there). Added argv[13] MODE: prior = Gaussian rho prior N(RHO,ERHO); free = uniform 0.1-60 rho_sun, no prior.
Data: flattened LC (prep3.py) within max(3*T14,0.25 d) of transit centre. Params t0,P,ln(Rp/R*),b,rho,q1,q2 (Kipping 2013, uniform),ln(jitter),f0; circular; batman supersampled x15 over 29.4 min.
ln L = -0.5*sum[(f-m)^2/s^2 + ln(2 pi s^2)], one white-noise s fitted; chi2_maxpost = sum((f-m)/s)^2 at max-posterior sample, N data, 9 free params.
usage: fit3.py EPIC SOURCE P T0 DUR_H RHO E_RHO RSTAR E_RSTAR TEFF E_TEFF NSTEP MODE"""
import sys, json, numpy as np, batman, emcee, time
from scipy.optimize import minimize
epic, src = sys.argv[1], sys.argv[2]
P0, T00, dur = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])/24
RHO, ERHO, RS, ERS, TEFF, ETEFF = map(float, sys.argv[6:12]); NSTEP = int(sys.argv[12]); MODE = sys.argv[13]
D = f'k2-deep/{epic}'
z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
ph = ((t - T00 + 0.5*P0) % P0) - 0.5*P0; w = np.abs(ph) < max(3*dur, 0.25)
t, f = t[w], f[w]
nmid = np.round(np.median(np.round((t - T00)/P0))); Tref = T00 + nmid*P0   # reference epoch near data middle
G = 6.674e-8; RSUN = 6.957e10; RHOSUN = 1.989e33/(4/3*np.pi*RSUN**3)
pm = batman.TransitParams(); pm.t0=0; pm.per=P0; pm.rp=0.02; pm.a=30; pm.inc=89.5; pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
m = batman.TransitModel(pm, t - Tref, supersample_factor=15, exp_time=1765.5/86400)
def a_rs(rho, P): return ((G*rho*RHOSUN*(P*86400)**2)/(3*np.pi))**(1/3)
def model(th):
    t0, P, lrp, b, rho, q1, q2, lsig, f0 = th
    a = a_rs(rho, P); inc = np.degrees(np.arccos(b/a))
    pm.t0 = t0; pm.per = P; pm.rp = np.exp(lrp); pm.a = a; pm.inc = inc
    pm.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
    return f0*m.light_curve(pm)
def lnp(th):
    t0, P, lrp, b, rho, q1, q2, lsig, f0 = th
    rp = np.exp(lrp)
    if not (-0.1<t0<0.1 and abs(P-P0)<0.05 and 0.001<rp<0.3 and 0<=b<1+rp and 0.1<rho<60 and 0<q1<1 and 0<q2<1 and -12<lsig<-5): return -np.inf
    if b/a_rs(rho,P) >= 1: return -np.inf
    lp = -0.5*((rho-RHO)/ERHO)**2 if MODE == "prior" else 0.0
    s2 = np.exp(2*lsig)
    r = f - model(th)
    return lp - 0.5*np.sum(r**2/s2 + np.log(2*np.pi*s2))
x0 = np.array([0., P0, np.log(0.023), 0.3, RHO, 0.4, 0.3, np.log(np.std(f)), 1.0])
nll = lambda x: -lnp(x)
best = minimize(nll, x0, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-6))
x1 = best.x if np.isfinite(best.fun) else x0
nw = 40; ndim = len(x1)
scl = np.array([1e-3, 1e-5, 0.02, 0.05, 0.3, 0.05, 0.05, 0.02, 1e-5])
p0 = []
while len(p0) < nw:
    c = x1 + scl*np.random.randn(ndim); c[3] = abs(c[3]); 
    if np.isfinite(lnp(c)): p0.append(c)
s = emcee.EnsembleSampler(nw, ndim, lnp)
t1 = time.time(); s.run_mcmc(np.array(p0), NSTEP, progress=False)
ch = s.get_chain(discard=NSTEP//3, thin=5, flat=True); lpch = s.get_log_prob(discard=NSTEP//3, thin=5, flat=True)
try: tau = s.get_autocorr_time(discard=NSTEP//3, quiet=True).tolist()
except Exception as e: tau = str(e)
t0, P, lrp, b, rho, q1, q2, lsig, f0 = ch.T
rp = np.exp(lrp); a = a_rs(rho, P); inc = np.arccos(b/a)
T14 = P/np.pi*np.arcsin(np.sqrt(np.clip((1+rp)**2 - b**2, 0, None))/(a*np.sin(inc)))*24
depth = []
idx = np.random.choice(len(ch), 2000, replace=False)
for i in idx:  # model depth at mid-transit, no supersampling
    pm2 = batman.TransitParams(); pm2.t0=0; pm2.per=P[i]; pm2.rp=rp[i]; pm2.a=a[i]; pm2.inc=np.degrees(inc[i]); pm2.ecc=0; pm2.w=90; pm2.limb_dark='quadratic'
    pm2.u=[2*np.sqrt(q1[i])*q2[i], np.sqrt(q1[i])*(1-2*q2[i])]
    depth.append(1 - batman.TransitModel(pm2, np.array([0.0])).light_curve(pm2)[0])
depth = np.array(depth)*1e6
Rs = RS + ERS*np.random.randn(len(ch)); Rp = rp*Rs*109.076
Te = TEFF + ETEFF*np.random.randn(len(ch)); Teq = Te*np.sqrt(1/(2*a))   # zero albedo, full redistribution
S = (Rs**2)*(Te/5772.)**4*(a*Rs*0.00465047)**-2
def q(x): p = np.percentile(x, [16, 50, 84]); return [float(p[1]), float(p[1]-p[0]), float(p[2]-p[1])]
res = dict(epic=epic, source=src, Tref_epoch=float(nmid), T0_BKJD=q(Tref + t0), P_d=q(P), RpRs=q(rp), b=q(b), rho_sun=q(rho), a_Rs=q(a),
           inc_deg=q(np.degrees(inc)), depth_ppm=q(depth), T14_h=q(T14), Rp_earth=q(Rp), Teq_K=q(Teq), S_earth=q(S), q1=q(q1), q2=q(q2),
           jitter_ppm=q(np.exp(lsig)*1e6), acc=float(np.mean(s.acceptance_fraction)), tau=tau, nstep=NSTEP, nwalk=nw, npts=int(len(t)),
           maxlnp=float(lpch.max()), priors=dict(rho=[RHO, ERHO], Rstar=[RS, ERS], Teff=[TEFF, ETEFF], ld='uniform q1,q2 (Kipping 2013)', ecc='fixed 0'),
           runtime_s=time.time()-t1, mode=MODE)
bi0 = ch[np.argmax(lpch)]; r0 = f - model(bi0); res['chi2_maxpost'] = float(np.sum((r0/np.exp(bi0[7]))**2)); res['ndata'] = int(len(f)); res['nfree'] = 9
json.dump(res, open(f'{D}/fit/result_{MODE}_{src}.json', 'w'), indent=1)
np.savez(f'{D}/fit/chain_{MODE}_{src}.npz', ch=ch, lp=lpch, Tref=Tref)
bi = ch[np.argmax(lpch)]; np.savez(f'{D}/fit/bestmodel_{MODE}_{src}.npz', t=t, f=f, mod=model(bi), Tref=Tref, best=bi)
for k, v in res.items(): print(k, v)
