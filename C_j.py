#!/usr/bin/env python3
"""
COMPLETE DEDUCTION OF C_J FOR EARTH-MOON TRANSFERS VIA L1
Based on analysis of Rust code (lib.rs, geocentric.rs, manifold.rs)

Author: astro_tfc v2.4.1 analysis
Date: June 2026
"""

import numpy as np

# ============================================================================
# 1. FUNDAMENTAL CONSTANTS OF THE EARTH-MOON SYSTEM
# ============================================================================
MU = 0.01215058560962404
D_CHAR = 384400000.0
V_CHAR = 1024.5
T_CHAR = D_CHAR / V_CHAR

print("=" * 80)
print("DEDUCTION OF C_J FOR EARTH-MOON TRANSFERS VIA L1")
print("=" * 80)
print("\nSTEP 1: SYSTEM CONSTANTS")
print(f"  mu = {MU:.15f}")
print(f"  D* = {D_CHAR:.0f} m")
print(f"  V* = {V_CHAR:.1f} m/s")
print(f"  T* = {T_CHAR:.0f} s = {T_CHAR/86400:.3f} days")

# ============================================================================
# 2. L1 POSITION
# ============================================================================
def calculate_L1_position(mu=MU, tol=1e-15):
    xi0 = (mu/3)**(1/3) * (1 - mu/3 - mu**2/9 - 23*mu**3/81)
    x = 1 - mu - xi0
    
    for _ in range(50):
        r1 = x + mu
        r2 = 1 - mu - x
        r1_3 = r1**3
        r2_3 = r2**3
        
        f = x - (1-mu)*r1/r1_3 - mu*(x-1+mu)/r2_3
        df = 1 + 2*(1-mu)/r1_3 + 2*mu/r2_3
        
        dx = f/df
        x -= dx
        if abs(dx) < tol:
            break
    
    return x, r1, r2

xL1, r1_L1, r2_L1 = calculate_L1_position()

print("\nSTEP 2: L1 POSITION")
print(f"  x_L1 = {xL1:.15f} (canonical)")
print(f"        = {xL1 * D_CHAR/1000:.0f} km from barycenter")
print(f"  r1 = {r1_L1:.6f} (to Earth, {r1_L1*D_CHAR/1000:.0f} km)")
print(f"  r2 = {r2_L1:.6f} (to Moon, {r2_L1*D_CHAR/1000:.0f} km)")

# ============================================================================
# 3. SECOND DERIVATIVES OF THE PSEUDO-POTENTIAL
# ============================================================================
def calculate_Omega_derivatives(xL1, mu=MU):
    r1 = xL1 + mu
    r2 = 1 - mu - xL1
    r1_3 = r1**3
    r2_3 = r2**3
    
    Omega_xx = 1 + 4*(1-mu)/r1_3 + 4*mu/r2_3
    Omega_yy = 1 - (1-mu)/r1_3 - mu/r2_3
    Omega_zz = -(1-mu)/r1_3 - mu/r2_3
    
    return Omega_xx, Omega_yy, Omega_zz

Omega_xx, Omega_yy, Omega_zz = calculate_Omega_derivatives(xL1)

print("\nSTEP 3: SECOND DERIVATIVES OF PSEUDO-POTENTIAL AT L1")
print(f"  Omega_xx = {Omega_xx:.6f} (confining)")
print(f"  Omega_yy = {Omega_yy:.6f} (repulsive)")
print(f"  Omega_zz = {Omega_zz:.6f} (repulsive)")

# ============================================================================
# 4. UNSTABLE EIGENVALUE
# ============================================================================
def calculate_eigenvalue(Omega_xx):
    a = 1 + 2*Omega_xx
    lambda_sq = (a + np.sqrt(a**2 - 4)) / 2
    lambda_u = np.sqrt(lambda_sq)
    return lambda_u, lambda_sq, a

lambda_u, lambda_sq, a = calculate_eigenvalue(Omega_xx)

print("\nSTEP 4: UNSTABLE EIGENVALUE OF L1")
print(f"  a = 1 + 2*Omega_xx = {a:.6f}")
print(f"  lambda^2 = {lambda_sq:.6f}")
print(f"  lambda = {lambda_u:.6f}")

# ============================================================================
# 5. EIGENVECTOR SHAPE FACTOR
# ============================================================================
def calculate_factor_k(lambda_u, Omega_yy):
    k = 2 * lambda_u / (lambda_u**2 - Omega_yy)
    return k

k = calculate_factor_k(lambda_u, Omega_yy)

print("\nSTEP 5: EIGENVECTOR SHAPE FACTOR")
print(f"  k = {k:.6f}")
print(f"  Unstable eigenvector: (1, {k:.3f}, 0, {lambda_u:.3f}, {lambda_u*k:.3f}, 0)")

# ============================================================================
# 6. SENSITIVITY alpha
# ============================================================================
def calculate_alpha(Omega_xx, Omega_yy, lambda_u, k):
    alpha = Omega_xx + k**2 * Omega_yy - lambda_u**2 * (1 + k**2)
    return alpha

alpha = calculate_alpha(Omega_xx, Omega_yy, lambda_u, k)

print("\nSTEP 6: SENSITIVITY alpha = dC_J/d(epsilon^2)")
print(f"  alpha = Omega_xx + k^2*Omega_yy - lambda^2*(1 + k^2)")
print(f"  alpha = {Omega_xx:.6f} + {k**2:.6f}*({Omega_yy:.6f}) - {lambda_u**2:.6f}*(1 + {k**2:.6f})")
print(f"  alpha = {alpha:.6f}")
print(f"  NOTE: alpha < 0 => C_J DECREASES with increasing epsilon")

# ============================================================================
# 7. JACOBI CONSTANT AT L1
# ============================================================================
def jacobi_constant(state, mu=MU):
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2) + (1-mu)/r1 + mu/r2
    v_sq = vx**2 + vy**2 + vz**2
    return 2*U - v_sq

CJ_L1 = jacobi_constant([xL1, 0, 0, 0, 0, 0])

print("\nSTEP 7: JACOBI CONSTANT AT L1")
print(f"  C_J(L1) = {CJ_L1:.10f}")
print(f"  To cross to the Moon: C_J < {CJ_L1:.6f}")

# ============================================================================
# 8. FUNDAMENTAL RELATIONSHIP C_J <-> epsilon
# ============================================================================
def epsilon_from_CJ(CJ_target):
    delta = CJ_L1 - CJ_target
    if delta <= 0:
        return 0.0
    return np.sqrt(delta / abs(alpha))

def CJ_from_epsilon(epsilon):
    return CJ_L1 + alpha * epsilon**2

print("\n" + "=" * 80)
print("STEP 8: FUNDAMENTAL RELATIONSHIP C_J <-> epsilon")
print("=" * 80)
print(f"\n  MASTER FORMULA:")
print(f"  C_J(epsilon) = C_J(L1) + alpha * epsilon^2")
print(f"  C_J(epsilon) = {CJ_L1:.6f} + ({alpha:.6f}) * epsilon^2")
print(f"\n  epsilon(C_J) = sqrt((C_J(L1) - C_J) / |alpha|)")
print(f"  epsilon(C_J) = sqrt(({CJ_L1:.6f} - C_J) / {abs(alpha):.6f})")

print("\n  TABLE C_J <-> epsilon:")
print(f"  {'C_J':<10} {'Delta C_J':<12} {'epsilon':<12} {'v_inj (m/s)':<15}")
print("  " + "-" * 49)

for cj in [3.188, 3.187, 3.186, 3.185, 3.184, 3.182, 3.180, 3.175, 3.170]:
    eps = epsilon_from_CJ(cj)
    delta = CJ_L1 - cj
    v_inj = eps * lambda_u * np.sqrt(1 + k**2) * V_CHAR
    print(f"  {cj:<10.4f} {delta:<12.6f} {eps:<12.6f} {v_inj:<15.1f}")

# ============================================================================
# 9. NEWTON-RAPHSON ALGORITHM
# ============================================================================
def simulate_newton_raphson(CJ_target, eps_guess=0.005, tol=1e-10, max_iter=30):
    direction = [-1.0, -k, 0.0, lambda_u, lambda_u * k, 0.0]
    eps = eps_guess
    history = []
    
    for i in range(max_iter):
        state = [
            xL1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        ]
        cj_test = jacobi_constant(state)
        err = cj_test - CJ_target
        
        history.append({
            'iter': i,
            'eps': eps,
            'CJ': cj_test,
            'err': err
        })
        
        if abs(err) < tol:
            return True, i+1, history
        
        eps_sq = eps * eps
        eps_sq_new = eps_sq - err / alpha
        
        if eps_sq_new <= 0.0:
            eps = eps * 0.5
            if eps < 1e-10:
                return False, i+1, history
            continue
        
        eps = np.clip(np.sqrt(eps_sq_new), 1e-8, 0.1)
    
    return False, max_iter, history

print("\n" + "=" * 80)
print("STEP 9: NEWTON-RAPHSON ALGORITHM")
print("=" * 80)
print("\n  Update formula:")
print("  epsilon^2_new = epsilon^2_current - (C_J(epsilon) - C_J_target) / alpha")
print("\n  SIMULATION:")

for cj_test in [3.188, 3.187, 3.185, 3.180]:
    print(f"\n  C_J = {cj_test:.3f} (Delta C_J = {CJ_L1 - cj_test:.6f}):")
    converged, iters, history = simulate_newton_raphson(cj_test)
    
    if converged:
        eps_final = history[-1]['eps']
        cj_final = history[-1]['CJ']
        dv = eps_final * lambda_u * np.sqrt(1 + k**2) * V_CHAR
        print(f"    CONVERGES in {iters} iterations")
        print(f"    epsilon = {eps_final:.6f}")
        print(f"    C_J = {cj_final:.10f}")
        print(f"    Estimated Delta V = {dv:.1f} m/s")
    else:
        print(f"    FAILS after {iters} iterations")

# ============================================================================
# 10. TOF AND DELTA V ESTIMATION
# ============================================================================
def estimate_TOF(epsilon):
    if epsilon < 1e-12:
        return float('inf')
    tau = (1.0 / lambda_u) * np.log(0.15 / epsilon)
    return tau * T_CHAR / 86400.0

def estimate_DV(epsilon):
    v_canonical = lambda_u * epsilon * np.sqrt(1 + k**2)
    return v_canonical * V_CHAR

print("\n" + "=" * 80)
print("STEP 10: TOF AND DELTA V ESTIMATION")
print("=" * 80)

# Empirical correction factors (calibrated with real data: 3.185 -> 40.4d, 67.1m/s)
eps_cal = epsilon_from_CJ(3.185)
tof_theoretical = estimate_TOF(eps_cal)
dv_theoretical = estimate_DV(eps_cal)
f_tof = 40.4 / tof_theoretical
f_dv = 67.1 / dv_theoretical

print(f"\n  Linear model of exponential escape:")
print(f"  TOF = (1/lambda) * ln(0.15/epsilon) * (T*/86400)")
print(f"  DV = lambda * epsilon * sqrt(1 + k^2) * V*")
print(f"\n  Calibration with real data (C_J=3.185):")
print(f"  Theoretical TOF = {tof_theoretical:.1f} days -> correction factor = {f_tof:.2f}")
print(f"  Theoretical DV = {dv_theoretical:.1f} m/s -> correction factor = {f_dv:.2f}")

print("\n  SWEEP OF FEASIBLE SOLUTIONS:")
print(f"  {'C_J':<10} {'epsilon':<12} {'TOF(d)':<12} {'DV(m/s)':<12} {'Product':<12}")
print("  " + "-" * 58)

best_product = float('inf')
best_cj = 0

for cj in np.arange(3.150, 3.1885, 0.0005):
    eps = epsilon_from_CJ(cj)
    if eps < 0.003:
        continue
    
    tof = estimate_TOF(eps) * f_tof
    dv = estimate_DV(eps) * f_dv
    product = tof * dv
    
    if tof < 60 and dv < 120 and product < best_product:
        best_product = product
        best_cj = cj
    
    if abs(cj - round(cj, 3)) < 0.0001:
        print(f"  {cj:<10.4f} {eps:<12.6f} {tof:<12.1f} {dv:<12.1f} {product:<12.0f}")

print(f"\n  DEDUCED OPTIMAL C_J: {best_cj:.4f}")
print(f"  Criterion: minimize TOF * DV")

# ============================================================================
# 11. FINAL SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("SUMMARY: HOW THE BACKEND DETERMINES C_J (lib.rs)")
print("=" * 80)
print(f"""
1. IT DOES NOT DETERMINE C_J AUTOMATICALLY
   The backend receives C_J as USER INPUT (config.target_jacobi)
   Default value: 3.187

2. VALIDATES THE RECEIVED C_J:
   - Must be < C_J(L1) = {CJ_L1:.6f}
   - Must be > C_J(L1) - 0.05 = {CJ_L1 - 0.05:.6f}
   - Valid window: C_J in [{CJ_L1 - 0.05:.4f}, {CJ_L1:.4f}]

3. CALCULATES epsilon VIA NEWTON-RAPHSON:
   epsilon^2_new = epsilon^2_current - (C_J(epsilon) - C_J_target) / alpha
   where alpha = {alpha:.6f}

4. CONSTRUCTS THE INJECTION STATE:
   x  = xL1 - epsilon = {xL1:.4f} - epsilon
   y  = -k * epsilon = -{k:.4f} * epsilon
   vx = lambda * epsilon = {lambda_u:.4f} * epsilon
   vy = lambda * k * epsilon = {lambda_u*k:.4f} * epsilon

5. VERIFIES FEASIBILITY:
   Stage 2: Propagates toward the Moon
   Success if altitude < 500 km in < 60 days

CONCLUSION:
  C_J is a DESIGN PARAMETER chosen by the user.
  The physics of the system (alpha, lambda, k) determines C_J <-> (TOF, DV).
  The optimal value arises from the tradeoff between time and fuel.
  
  For your mission: C_J = 3.185 is an EXCELLENT choice.
  Theoretical optimal C_J: {best_cj:.4f}
  Difference: {abs(3.185 - best_cj):.4f}
""")

print("=" * 80)
print("ANALYSIS COMPLETED")
print("=" * 80)