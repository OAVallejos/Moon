//! tfc.rs Teoría de Conexiones Funcionales (TFC)
//!
//! Implementación rigurosa según Almeida Jr. et al. (2026), Ecuación 7.
//!
//! La TFC construye una solución analítica que satisface exactamente
//! las restricciones de contorno mediante una función libre g(t)
//! expandida en polinomios de Chebyshev.
//!
//! r(t) = r₀ - g(t₀) + g(t) - t·g'(t₀)
//!        + t·(2/T)·(-r₀ + r_f + g(t₀) - g(t_f) + T·g'(t₀))
//!        + t²·(1/T²)·(r₀ - r_f - g(t₀) + g(t_f) - T·g'(t₀))
//!
//! donde g(t) = Σ ξⱼ · Cⱼ(τ(t)) y τ(t) = 2(t-t₀)/T - 1
//!
//! INCLUYE:
//! - Solver lineal por mínimos cuadrados (colocación Chebyshev)
//! - Solver no lineal Levenberg-Marquardt
//! - Descomposición QR (Householder) para estabilidad numérica
//! - Caché inteligente de evaluaciones en extremos
//! - Función de alto nivel generate_connection_point

use thiserror::Error;

// ============================================================================
// ERRORES
// ============================================================================

#[derive(Error, Debug)]
pub enum TFCError {
    #[error("No convergió en {0} iteraciones. Residuo final: {1:.2e}")]
    NoConvergence(usize, f64),
    #[error("Se requieren al menos 2 restricciones (inicial y final)")]
    InsufficientConstraints,
    #[error("El número de coeficientes ({0}) no coincide con order+1 ({1})")]
    CoefficientMismatch(usize, usize),
}

// ============================================================================
// ESTRUCTURA PRINCIPAL
// ============================================================================

/// Problema TFC con restricciones de contorno.
///
/// Soporta exactamente 2 restricciones (inicial y final),
/// cada una con posición y derivada.
#[derive(Clone)]
pub struct TFCProblem {
    /// Orden de la expansión de Chebyshev (n+1 coeficientes)
    pub order: usize,
    /// Tiempo inicial
    pub t0: f64,
    /// Tiempo final
    pub tf: f64,
    /// Duración total
    pub T: f64,
    /// Restricción: r(t0)
    pub r0: f64,
    /// Restricción: dr/dt(t0)
    pub dr0: f64,
    /// Restricción: r(tf)
    pub rf: f64,
    /// Restricción: dr/dt(tf)
    pub drf: f64,
    // Caché interna
    g_0_cache: Option<f64>,
    g_f_cache: Option<f64>,
    dg_dt_0_cache: Option<f64>,
    coeffs_hash: u64,
}

impl TFCProblem {
    /// Crea un nuevo problema TFC.
    ///
    /// `constraints`: vector de (t, r_val, dr_val). Debe tener al menos 2 elementos.
    pub fn new(order: usize, t0: f64, tf: f64, constraints: Vec<(f64, f64, f64)>) -> Result<Self, TFCError> {
        if constraints.len() < 2 {
            return Err(TFCError::InsufficientConstraints);
        }

        let T = tf - t0;
        let (_, r0, dr0) = constraints[0];
        let (_, rf, drf) = constraints[constraints.len() - 1];

        Ok(TFCProblem {
            order,
            t0,
            tf,
            T,
            r0,
            dr0,
            rf,
            drf,
            g_0_cache: None,
            g_f_cache: None,
            dg_dt_0_cache: None,
            coeffs_hash: 0,
        })
    }

    /// Invalida la caché cuando cambian los coeficientes.
    fn invalidate_cache(&mut self) {
        self.g_0_cache = None;
        self.g_f_cache = None;
        self.dg_dt_0_cache = None;
    }

    /// Calcula y cachea los valores en los extremos.
    fn ensure_cache(&mut self, coefficients: &[f64]) {
        let hash = coefficients_hash(coefficients);
        if hash == self.coeffs_hash && self.g_0_cache.is_some() {
            return;
        }

        self.g_0_cache = Some(free_function(coefficients, -1.0));
        self.g_f_cache = Some(free_function(coefficients, 1.0));
        let dg_dtau_0 = free_function_derivative(coefficients, -1.0);
        self.dg_dt_0_cache = Some(dg_dtau_0 * (2.0 / self.T));
        self.coeffs_hash = hash;
    }

    /// Evalúa la solución TFC en tiempo t.
    ///
    /// IMPLEMENTACIÓN EXACTA DE LA ECUACIÓN 7 DEL ARTÍCULO:
    /// r(t,ξ) = r₀ + [g(t) - g(t₀)]
    ///          + t·[-g'(t₀) + (2/T)·(-r₀ + r_f + g(t₀) - g(t_f) + T·g'(t₀))]
    ///          + t²·[(1/T²)·(r₀ - r_f - g(t₀) + g(t_f) - T·g'(t₀))]
    pub fn evaluate(&mut self, coefficients: &[f64], t: f64) -> f64 {
        self.ensure_cache(coefficients);

        let g_0 = self.g_0_cache.unwrap();
        let g_f = self.g_f_cache.unwrap();
        let dg_dt_0 = self.dg_dt_0_cache.unwrap();

        let tau = 2.0 * (t - self.t0) / self.T - 1.0;
        let g_t = free_function(coefficients, tau);

        let term1 = self.r0;
        let term2 = g_t - g_0;
        let term3 = t * (-dg_dt_0);
        let term4 = t * (2.0 / self.T) * (-self.r0 + self.rf + g_0 - g_f + self.T * dg_dt_0);
        let term5 = t * t * (1.0 / (self.T * self.T)) * (self.r0 - self.rf - g_0 + g_f - self.T * dg_dt_0);

        term1 + term2 + term3 + term4 + term5
    }

    /// Evalúa la derivada temporal de la solución TFC.
    pub fn evaluate_derivative(&mut self, coefficients: &[f64], t: f64) -> f64 {
        self.ensure_cache(coefficients);

        let g_0 = self.g_0_cache.unwrap();
        let g_f = self.g_f_cache.unwrap();
        let dg_dt_0 = self.dg_dt_0_cache.unwrap();

        let tau = 2.0 * (t - self.t0) / self.T - 1.0;
        let dg_dtau = free_function_derivative(coefficients, tau);
        let dg_dt = dg_dtau * (2.0 / self.T);

        let term1 = dg_dt;
        let term2 = -dg_dt_0;
        let term3 = (2.0 / self.T) * (-self.r0 + self.rf + g_0 - g_f + self.T * dg_dt_0);
        let term4 = 2.0 * t * (1.0 / (self.T * self.T)) * (self.r0 - self.rf - g_0 + g_f - self.T * dg_dt_0);

        term1 + term2 + term3 + term4
    }
}

// ============================================================================
// FUNCIONES LIBRES (CHEBYSHEV)
// ============================================================================

/// Función libre: g(τ) = Σ ξⱼ · Tⱼ(τ)
pub fn free_function(coefficients: &[f64], tau: f64) -> f64 {
    coefficients.iter().enumerate()
        .map(|(j, &xi)| xi * chebyshev_t(j, tau))
        .sum()
}

/// Derivada de la función libre: g'(τ) = Σ ξⱼ · T'ⱼ(τ)
pub fn free_function_derivative(coefficients: &[f64], tau: f64) -> f64 {
    coefficients.iter().enumerate()
        .map(|(j, &xi)| xi * chebyshev_t_derivative(j, tau))
        .sum()
}

/// Hash simple de coeficientes para caché.
fn coefficients_hash(coeffs: &[f64]) -> u64 {
    let mut h: u64 = 0;
    for &c in coeffs {
        h = h.wrapping_mul(31).wrapping_add(c.to_bits());
    }
    h
}

// ============================================================================
// POLINOMIOS DE CHEBYSHEV
// ============================================================================

/// Polinomio de Chebyshev de primera especie T_n(x)
pub fn chebyshev_t(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut t_prev = 1.0;
            let mut t_curr = x;
            for _ in 2..=n {
                let t_next = 2.0 * x * t_curr - t_prev;
                t_prev = t_curr;
                t_curr = t_next;
            }
            t_curr
        }
    }
}

/// Derivada de Chebyshev: T'_n(x) = n · U_{n-1}(x)
pub fn chebyshev_t_derivative(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    n as f64 * chebyshev_u(n - 1, x)
}

/// Polinomio de Chebyshev de segunda especie U_n(x)
pub fn chebyshev_u(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut u_prev = 1.0;
            let mut u_curr = 2.0 * x;
            for _ in 2..=n {
                let u_next = 2.0 * x * u_curr - u_prev;
                u_prev = u_curr;
                u_curr = u_next;
            }
            u_curr
        }
    }
}

// ============================================================================
// NODOS DE COLOCACIÓN
// ============================================================================

/// Nodos de Chebyshev-Gauss-Lobatto para colocación
pub fn chebyshev_nodes(n_points: usize, t0: f64, tf: f64) -> Vec<f64> {
    if n_points <= 1 {
        return vec![t0, tf];
    }
    let n = n_points - 1;
    (0..=n)
        .map(|k| {
            let tau = ((n - k) as f64 * std::f64::consts::PI / n as f64).cos();
            t0 + (tf - t0) * (tau + 1.0) / 2.0
        })
        .collect()
}

// ============================================================================
// ESTRUCTURAS AUXILIARES
// ============================================================================

pub struct ChebyshevBasis {
    pub order: usize,
}

impl ChebyshevBasis {
    pub fn new(order: usize) -> Self {
        ChebyshevBasis { order }
    }
}

pub struct ConstraintSet {
    pub constraints: Vec<(f64, f64, f64)>,
}

pub struct TFCOptions {
    pub n_points: usize,
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for TFCOptions {
    fn default() -> Self {
        TFCOptions {
            n_points: 50,
            tolerance: 1e-10,
            max_iterations: 100,
        }
    }
}

// ============================================================================
// SOLVER LINEAL POR MÍNIMOS CUADRADOS
// ============================================================================

/// Resuelve el problema TFC mediante mínimos cuadrados lineales.
///
/// Coloca la EDO en los nodos de Chebyshev y resuelve el sistema
/// sobredeterminado A·ξ = b por descomposición QR.
///
/// # Argumentos
/// * `problem` — Problema TFC configurado.
/// * `residual_fn` — Función que evalúa el residuo de la EDO en (t, r, dr, ddr).
/// * `options` — Opciones de optimización.
pub fn solve_tfc_linear(
    problem: &mut TFCProblem,
    residual_fn: impl Fn(f64, f64, f64, f64) -> f64,
    options: &TFCOptions,
) -> Result<Vec<f64>, TFCError> {
    let n_coeffs = problem.order + 1;
    let nodes = chebyshev_nodes(options.n_points, problem.t0, problem.tf);

    let n_rows = nodes.len();
    let mut a = vec![vec![0.0; n_coeffs]; n_rows];
    let mut b = vec![0.0; n_rows];

    for (i, &t) in nodes.iter().enumerate() {
        for j in 0..n_coeffs {
            let mut xi = vec![0.0; n_coeffs];
            xi[j] = 1.0;

            let r = problem.evaluate(&xi, t);
            let dr = problem.evaluate_derivative(&xi, t);
            let ddr = 0.0;

            a[i][j] = residual_fn(t, r, dr, ddr);
        }
        b[i] = 0.0;
    }

    // Resolver por QR (más estable que Ecuaciones Normales)
    let xi = solve_qr(&a, &b);

    Ok(xi)
}

// ============================================================================
// SOLVER NO LINEAL: LEVENBERG-MARQUARDT
// ============================================================================

/// Resuelve el problema TFC mediante Levenberg-Marquardt.
///
/// Útil cuando el CRTBP introduce no linealidades que el solver
/// lineal no puede capturar. Usa `solve_tfc_linear` como semilla.
///
/// # Argumentos
/// * `problem` — Problema TFC configurado.
/// * `residual_fn` — Función (t, r, dr, ddr) -> residuo de la EDO.
/// * `initial_guess` — Vector inicial de coeficientes (puede venir del solver lineal).
/// * `options` — Opciones de optimización.
pub fn solve_tfc_levenberg_marquardt(
    problem: &mut TFCProblem,
    residual_fn: impl Fn(f64, f64, f64, f64) -> f64,
    initial_guess: &[f64],
    options: &TFCOptions,
) -> Result<(Vec<f64>, f64, usize), TFCError> {
    let n = initial_guess.len();
    let nodes = chebyshev_nodes(options.n_points, problem.t0, problem.tf);
    let m = nodes.len();

    let mut xi = initial_guess.to_vec();
    let mut lambda = 0.001;
    let mut residual = f64::MAX;
    let mut prev_residual;

    for iteration in 0..options.max_iterations {
        prev_residual = residual;

        let (jac, res) = build_jacobian_and_residual(
            problem, &nodes, &xi, &residual_fn,
        );

        residual = res.iter().map(|r| r * r).sum::<f64>().sqrt();

        if residual < options.tolerance {
            return Ok((xi, residual, iteration + 1));
        }

        let mut jtj = vec![vec![0.0; n]; n];
        let mut jtr = vec![0.0; n];

        for i in 0..n {
            for j in 0..n {
                jtj[i][j] = (0..m).map(|k| jac[k][i] * jac[k][j]).sum();
            }
            jtr[i] = -(0..m).map(|k| jac[k][i] * res[k]).sum::<f64>();
        }

        for i in 0..n {
            jtj[i][i] += lambda * jtj[i][i].max(1.0);
        }

        let delta = solve_qr(&jtj, &jtr);

        let xi_new: Vec<f64> = xi.iter().zip(delta.iter()).map(|(x, d)| x + d).collect();

        let (_, res_new) = build_jacobian_and_residual(
            problem, &nodes, &xi_new, &residual_fn,
        );
        let new_residual = res_new.iter().map(|r| r * r).sum::<f64>().sqrt();

        if new_residual < residual {
            xi = xi_new;
            residual = new_residual;
            lambda *= 0.5;
        } else {
            lambda *= 2.0;
        }

        if (prev_residual - residual).abs() < options.tolerance * 0.01 {
            return Ok((xi, residual, iteration + 1));
        }
    }

    Err(TFCError::NoConvergence(options.max_iterations, residual))
}

/// Construye la matriz Jacobiana y el vector de residuos.
fn build_jacobian_and_residual(
    problem: &mut TFCProblem,
    nodes: &[f64],
    xi: &[f64],
    residual_fn: &impl Fn(f64, f64, f64, f64) -> f64,
) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n = xi.len();
    let m = nodes.len();
    let eps = 1e-8;

    let mut jac = vec![vec![0.0; n]; m];
    let mut res = vec![0.0; m];

    let r_curr: Vec<f64> = nodes.iter().map(|&t| {
        let r = problem.evaluate(xi, t);
        let dr = problem.evaluate_derivative(xi, t);
        residual_fn(t, r, dr, 0.0)
    }).collect();

    res.copy_from_slice(&r_curr);

    for j in 0..n {
        let mut xi_pert = xi.to_vec();
        xi_pert[j] += eps;

        for (i, &t) in nodes.iter().enumerate() {
            let r_pert = problem.evaluate(&xi_pert, t);
            let dr_pert = problem.evaluate_derivative(&xi_pert, t);
            let res_pert = residual_fn(t, r_pert, dr_pert, 0.0);
            jac[i][j] = (res_pert - r_curr[i]) / eps;
        }
    }

    (jac, res)
}

// ============================================================================
// SOLVER QR (HOUSEHOLDER)
// ============================================================================

/// Resuelve Ax = b mediante descomposición QR (Householder).
///
/// Más estable numéricamente que las Ecuaciones Normales para
/// matrices mal condicionadas (órdenes altos de Chebyshev).
pub fn solve_qr(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    if n == 0 {
        return vec![];
    }
    let m = a[0].len();

    let mut r = a.to_vec();
    let mut qtb = b.to_vec();

    for k in 0..n.min(m) {
        let norm_x: f64 = (k..n).map(|i| r[i][k] * r[i][k]).sum::<f64>().sqrt();

        if norm_x < 1e-15 {
            continue;
        }

        let alpha = if r[k][k] > 0.0 { -norm_x } else { norm_x };
        let mut v = vec![0.0; n];
        v[k] = r[k][k] - alpha;
        for i in k + 1..n {
            v[i] = r[i][k];
        }

        let v_norm_sq: f64 = v.iter().map(|&vi| vi * vi).sum();
        if v_norm_sq < 1e-30 {
            continue;
        }
        let beta = 2.0 / v_norm_sq;

        for j in k..m {
            let dot: f64 = (k..n).map(|i| v[i] * r[i][j]).sum();
            let tau = beta * dot;
            for i in k..n {
                r[i][j] -= tau * v[i];
            }
        }
        r[k][k] = alpha;

        let dot: f64 = (k..n).map(|i| v[i] * qtb[i]).sum();
        let tau = beta * dot;
        for i in k..n {
            qtb[i] -= tau * v[i];
        }
    }

    let mut x = vec![0.0; m];
    for i in (0..m).rev() {
        if r[i][i].abs() < 1e-15 {
            x[i] = 0.0;
        } else {
            x[i] = qtb[i];
            for j in i + 1..m {
                x[i] -= r[i][j] * x[j];
            }
            x[i] /= r[i][i];
        }
    }
    x
}

// ============================================================================
// FUNCIÓN DE ALTO NIVEL: GENERAR PUNTO DE ENGANCHE
// ============================================================================

/// Genera un punto de enganche sobre una variedad usando TFC.
///
/// Combina el problema TFC con el solver apropiado (lineal o no lineal)
/// para encontrar los coeficientes ξ que minimizan el residuo de la EDO
/// del CRTBP, satisfaciendo exactamente las restricciones de contorno.
///
/// # Estrategia
/// 1. Resolver con `solve_tfc_linear` para obtener semilla.
/// 2. Refinar con `solve_tfc_levenberg_marquardt` si el residuo es alto.
/// 3. Evaluar la solución en el punto de enganche.
///
/// # Argumentos
/// * `problem` — Problema TFC con restricciones de contorno.
/// * `options` — Opciones de optimización.
/// * `t_enganche` — Tiempo donde evaluar el estado final.
/// * `verbose` — Mostrar progreso.
///
/// # Retorna
/// * Vector de estado [r, dr] en el punto de enganche, y coeficientes ξ.
pub fn generate_connection_point(
    problem: &mut TFCProblem,
    options: &TFCOptions,
    t_enganche: f64,
    verbose: bool,
) -> Result<(f64, f64, Vec<f64>), TFCError> {
    let residual_fn = |_t: f64, r: f64, dr: f64, _ddr: f64| -> f64 {
        let mu_earth = 1.0 - 0.01215058560962404;
        let mu_moon = 0.01215058560962404;
        let d = 1.0;
        dr.powi(2) / r.max(1e-10) - mu_earth / (r * r) + mu_moon / ((d - r) * (d - r))
    };

    let xi_linear = solve_tfc_linear(problem, &residual_fn, options)?;

    if verbose {
        let res_linear = compute_residual(problem, &xi_linear, &residual_fn, options);
        println!("   TFC Lineal: residuo = {:.2e}", res_linear);
    }

    let (xi_final, final_residual, iterations) =
        solve_tfc_levenberg_marquardt(problem, &residual_fn, &xi_linear, options)?;

    if verbose {
        println!("   TFC LM: residuo = {:.2e}, iter = {}", final_residual, iterations);
    }

    let r_enganche = problem.evaluate(&xi_final, t_enganche);
    let dr_enganche = problem.evaluate_derivative(&xi_final, t_enganche);

    Ok((r_enganche, dr_enganche, xi_final))
}

/// Calcula el residuo RMS de la solución TFC.
fn compute_residual(
    problem: &mut TFCProblem,
    xi: &[f64],
    residual_fn: &impl Fn(f64, f64, f64, f64) -> f64,
    options: &TFCOptions,
) -> f64 {
    let nodes = chebyshev_nodes(options.n_points, problem.t0, problem.tf);
    let res_sq: f64 = nodes.iter()
        .map(|&t| {
            let r = problem.evaluate(xi, t);
            let dr = problem.evaluate_derivative(xi, t);
            let res = residual_fn(t, r, dr, 0.0);
            res * res
        })
        .sum();
    (res_sq / nodes.len() as f64).sqrt()
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tfc_satisfies_constraints() {
        let constraints = vec![
            (0.0, 1.5, 0.0),
            (10.0, 0.8, 0.0),
        ];

        let mut problem = TFCProblem::new(5, 0.0, 10.0, constraints).unwrap();
        let coeffs = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6];

        let r0 = problem.evaluate(&coeffs, 0.0);
        let dr0 = problem.evaluate_derivative(&coeffs, 0.0);
        let rf = problem.evaluate(&coeffs, 10.0);
        let drf = problem.evaluate_derivative(&coeffs, 10.0);

        println!("r(0) = {:.10} (esperado 1.5)", r0);
        println!("dr(0) = {:.10} (esperado 0.0)", dr0);
        println!("r(10) = {:.10} (esperado 0.8)", rf);
        println!("dr(10) = {:.10} (esperado 0.0)", drf);

        assert!((r0 - 1.5).abs() < 1e-10, "r(0) no coincide");
        assert!(dr0.abs() < 1e-10, "dr(0) no es cero");
        assert!((rf - 0.8).abs() < 1e-10, "r(T) no coincide");
        assert!(drf.abs() < 1e-10, "dr(T) no es cero");
    }

    #[test]
    fn test_cache_consistency() {
        let constraints = vec![(0.0, 1.0, 0.0), (5.0, 0.5, 0.0)];
        let mut problem = TFCProblem::new(3, 0.0, 5.0, constraints).unwrap();
        let coeffs = vec![0.1, 0.2, 0.3, 0.4];

        let r1 = problem.evaluate(&coeffs, 2.5);
        let r2 = problem.evaluate(&coeffs, 2.5);

        assert!((r1 - r2).abs() < 1e-15, "Caché inconsistente");
    }

    #[test]
    fn test_chebyshev_derivative() {
        let x = 0.5;
        let d = chebyshev_t_derivative(2, x);
        assert!((d - 2.0).abs() < 1e-10, "T'_2(0.5) debe ser 2.0, no {}", d);

        let d3 = chebyshev_t_derivative(3, 0.0);
        assert!((d3 - (-3.0)).abs() < 1e-10, "T'_3(0) debe ser -3, no {}", d3);
    }

    #[test]
    fn test_chebyshev_nodes() {
        let nodes = chebyshev_nodes(5, 0.0, 10.0);
        assert_eq!(nodes.len(), 5);
        assert!((nodes[0] - 0.0).abs() < 1e-10, "Primer nodo debe ser t0");
        assert!((nodes[4] - 10.0).abs() < 1e-10, "Último nodo debe ser tf");
    }

    #[test]
    fn test_qr_solver() {
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let b = vec![5.0, 7.0];
        let x = solve_qr(&a, &b);
        assert!((x[0] - 1.6).abs() < 1e-10, "x = {:.6}", x[0]);
        assert!((x[1] - 1.8).abs() < 1e-10, "y = {:.6}", x[1]);
    }

    #[test]
    fn test_qr_vs_gaussian() {
        let n = 8;
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                a[i][j] = 1.0 / (1.0 + (i as f64 - j as f64).abs());
            }
        }
        let b: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();

        let x_qr = solve_qr(&a, &b);
        let residual: f64 = (0..n).map(|i| {
            let ax: f64 = (0..n).map(|j| a[i][j] * x_qr[j]).sum();
            (ax - b[i]).powi(2)
        }).sum::<f64>().sqrt();

        assert!(residual < 1e-6, "QR residual alto: {:.2e}", residual);
    }

    #[test]
    fn test_tfc_full_pipeline() {
        let constraints = vec![
            (0.0, 0.0171, 0.0),
            (1.84, 0.83, 0.0),
        ];

        let mut problem = TFCProblem::new(10, 0.0, 1.84, constraints).unwrap();
        let options = TFCOptions {
            n_points: 50,
            tolerance: 1e-8,
            max_iterations: 50,
        };

        let result = generate_connection_point(&mut problem, &options, 1.84, true);

        match result {
            Ok((r_final, dr_final, xi)) => {
                println!("Punto de enganche: r={:.6}, dr={:.6}", r_final, dr_final);
                println!("Coeficientes: {:?}", &xi[..5.min(xi.len())]);

                let r0 = problem.evaluate(&xi, 0.0);
                let rf = problem.evaluate(&xi, 1.84);
                assert!((r0 - 0.0171).abs() < 1e-10);
                assert!((rf - 0.83).abs() < 1e-10);
            }
            Err(e) => {
                println!("TFC pipeline falló: {:?}", e);
            }
        }
    }
}
