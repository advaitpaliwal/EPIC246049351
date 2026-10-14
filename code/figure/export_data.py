# Export the detrended EVEREST light curve (with phase and median model) and a posterior summary as CSV. usage: python -I export_data.py K2DEEP_DIR OUTDIR
import sys, json, numpy as np, batman
D, O = sys.argv[1], sys.argv[2]
R = json.load(open(f'{D}/fit/result_prior_everest.json')); z = np.load(f'{D}/lc/everest.npz')
P, T0 = R['P_d'][0], R['T0_BKJD'][0]; q1, q2 = R['q1'][0], R['q2'][0]
p = batman.TransitParams(); p.t0 = T0; p.per = P; p.rp = R['RpRs'][0]; p.a = R['a_Rs'][0]; p.inc = R['inc_deg'][0]
p.ecc = 0.; p.w = 90.; p.limb_dark = 'quadratic'; p.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
t, f = z['t'], z['f']; m = batman.TransitModel(p, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(p)
ph = (t - T0 + 0.5*P) % P - 0.5*P
np.savetxt(f'{O}/lc/everest_detrended_c12.csv', np.column_stack([t, t + 2454833, ph*24, (f - 1)*1e6, (m - 1)*1e6]), delimiter=',', fmt='%.7f',
           header='time_bkjd,time_bjd_tdb,hours_from_midtransit,flux_ppm,median_model_ppm', comments='')
np.savetxt(f'{O}/lc/everest_raw_c12.csv', np.column_stack([z['t_raw'], z['f_raw']]), delimiter=',', fmt='%.7f', header='time_bkjd,flux_normalised', comments='')
keys = ['T0_BKJD', 'P_d', 'RpRs', 'b', 'rho_sun', 'a_Rs', 'inc_deg', 'depth_ppm', 'T14_h', 'Rp_earth', 'Teq_K', 'S_earth', 'q1', 'q2', 'jitter_ppm']
with open(f'{O}/fit/posterior_summary_everest.csv', 'w') as fh:
    fh.write('parameter,median,minus_1sigma,plus_1sigma\n')
    for k in keys: fh.write('%s,%.8g,%.4g,%.4g\n' % (k, *R[k]))
print('ok', len(t))
