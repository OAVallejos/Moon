# test_permutacion_vcv_filtrado.py
#!/usr/bin/env python3     
import numpy as np         import glob, os

print("=" * 80)            
print("  PERMUTATION TEST WITH DYNAMIC UV FILTER (ONLY 50% LONGEST BASELINES)")
print("=" * 80)

archivos_npz = sorted(glob.glob('output/stokes/*_vis.npz'))
num_bins = 360

for arch in archivos_npz:
    nombre_base = os.path.basename(arch).replace('_vis.npz', '')

    with np.load(arch) as data:
        u, v = data['u'], data['v']
        vis_q, vis_u = data['vis_q'], data['vis_u']

    uv_dist = np.sqrt(u**2 + v**2)

    # DYNAMIC FILTER: We cut by the median to isolate high spatial frequencies
    UV_CORTE = np.median(uv_dist)
    mask_uv = uv_dist >= UV_CORTE

    u_filt = u[mask_uv]
    v_filt = v[mask_uv]
    vis_q_filt = vis_q[mask_uv]
    vis_u_filt = vis_u[mask_uv]

    def calcular_espectro(u_arr, v_arr, q_arr, u_arr_vis, barajar=False):
        if barajar:
            # Independent permutation to break Stokes pairs at the root
            idx_q = np.random.permutation(len(u_arr))
            idx_u = np.random.permutation(len(u_arr))
            q_arr = q_arr[idx_q]
            u_arr_vis = u_arr_vis[idx_u]

        phi = np.arctan2(v_arr, u_arr)
        phi = np.where(phi < 0, phi + 2*np.pi, phi)
        chi = 0.5 * np.arctan2(u_arr_vis, q_arr)

        bins = np.linspace(0, 2*np.pi, num_bins+1)
        perfil = np.zeros(num_bins)
        counts = np.zeros(num_bins)
        idx_bin = np.digitize(phi, bins) - 1

        for i in range(len(phi)):
            if 0 <= idx_bin[i] < num_bins:
                perfil[idx_bin[i]] += chi[i]
                counts[idx_bin[i]] += 1

        with np.errstate(divide='ignore', invalid='ignore'):
            perfil_prom = np.where(counts > 5, perfil/counts, 0.0)

        fft = np.fft.fft(perfil_prom)
        power = np.abs(fft)**2 / num_bins
        return power[2], power[6]

    m2_orig, m6_orig = calcular_espectro(u_filt, v_filt, vis_q_filt, vis_u_filt, False)

    np.random.seed(42)
    m2_bar, m6_bar = calcular_espectro(u_filt, v_filt, vis_q_filt, vis_u_filt, True)

    frac_m6 = m6_bar / m6_orig if m6_orig > 1e-10 else 0.0

    print(f"\n{nombre_base}")
    print(f"  Points analyzed (UV cutoff at {UV_CORTE:.4f} light-s): {len(u_filt)}/{len(u)}")
    print(f"  m=2 [Orig]: {m2_orig:.4e} | [Shuffled]: {m2_bar:.4e}")
    print(f"  m=6 [Orig]: {m6_orig:.4e} | [Shuffled]: {m6_bar:.4e}")

    if frac_m6 < 0.5:
        print(f"  ✅ m=6 COLLAPSES on shuffling (fraction: {frac_m6:.4f}) → PHYSICAL SIGNAL")
    else:
        print(f"  ❌ m=6 SURVIVES shuffling (fraction: {frac_m6:.4f}) → INSTRUMENTAL CONTAMINATION")
    print("-" * 80)
    