# test_permutation_vcv.py
#!/usr/bin/env python3     
import numpy as np         import glob                
import os                                             print("======================================================================")  print("     STOCHASTIC PERMUTATION TEST (SHUFFLE) - MODES M2 AND M6")
print("======================================================================\n")

archivos_npz = sorted(glob.glob('output/stokes/*_vis.npz'))
num_bins = 360
np.random.seed(42)  # Fixed seed for scientific reproducibility

for arch in archivos_npz:
    nombre_base = os.path.basename(arch).replace('_vis.npz', '')

    with np.load(arch) as data:
        u, v = data['u'], data['v']
        vis_q, vis_u = data['vis_q'], data['vis_u']

    phi_orig = np.arctan2(v, u)
    phi_orig = np.where(phi_orig < 0, phi_orig + 2 * np.pi, phi_orig)
    chi_uv = 0.5 * np.arctan2(vis_u, vis_q)

    # THE TRUE DECOUPLING: We randomly shuffle the EVPA
    # This destroys the spatial coherence of the source
    chi_barajado = np.random.permutation(chi_uv)

    bin_edges = np.linspace(0, 2 * np.pi, num_bins + 1)
    idx = np.digitize(phi_orig, bin_edges) - 1

    p_orig, c_orig = np.zeros(num_bins), np.zeros(num_bins)
    p_shuf, c_shuf = np.zeros(num_bins), np.zeros(num_bins)

    for i in range(len(phi_orig)):
        if 0 <= idx[i] < num_bins:
            p_orig[idx[i]] += chi_uv[i]
            c_orig[idx[i]] += 1.0

            p_shuf[idx[i]] += chi_barajado[i]
            c_shuf[idx[i]] += 1.0

    perfil_orig = np.zeros(num_bins)
    perfil_shuf = np.zeros(num_bins)
    np.divide(p_orig, c_orig, out=perfil_orig, where=c_orig > 0)
    np.divide(p_shuf, c_shuf, out=perfil_shuf, where=c_shuf > 0)

    esp_orig = np.abs(np.fft.fft(perfil_orig))**2 / num_bins
    esp_shuf = np.abs(np.fft.fft(perfil_shuf))**2 / num_bins

    print(f"Dataset: {nombre_base}")
    print(f"  -> m=2 [Original]: {esp_orig[2]:.4e} | [Shuffled]: {esp_shuf[2]:.4e}")
    print(f"  -> m=6 [Original]: {esp_orig[6]:.4e} | [Shuffled]: {esp_shuf[6]:.4e}")

    ratio_m6 = esp_shuf[6] / esp_orig[6] if esp_orig[6] > 0 else 0
    print(f"  -> Remaining m=6 power fraction: {ratio_m6:.4f}")
    print("-" * 75)