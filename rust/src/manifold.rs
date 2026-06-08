//! Variedades invariantes y punto L1
//! Basado en Almeida Jr. et al. (2026) - Tabla 4
//! 
//! Implementación estrictamente analítica con diagonalización exacta de la matriz
//! jacobiana en L1. Cero valores hardcodeados. Cero diferencias finitas.
//! 
//! Los autovectores se derivan de las derivadas segundas del pseudo-potencial
//! efectivo Ω. La sensibilidad de la constante de Jacobi se calcula mediante
//! la forma cuadrática exacta del Hessiano, sin perturbaciones numéricas.
//!
//! # Constantes geométricas precalculadas
//!
//! La sensibilidad α = dC_J/d(ε²) evaluada en L1 es una constante del sistema
//! Tierra-Luna. No depende de la dirección (el cuadrado del autovector elimina
//! el signo) ni del target_cj. Se expone mediante `jacobi_sensitivity_l1()`.

use thiserror::Error;
use crate::crtbp::{StateVector, jacobi_constant, MU};

// ============================================================================
// PUNTO DE EQUILIBRIO L1
// ============================================================================

/// Resuelve la posición exacta de L1 mediante Newton-Raphson.
///
/// Ecuación de equilibrio: dΩ/dx = 0 en el segmento Tierra-Luna.
///
/// El pseudo-potencial efectivo en el sistema rotante es:
///
/// Ω(x, y, z) = (1/2)(x² + y²) + (1-μ)/r_e + μ/r_m
///
/// donde r_e = distancia a la Tierra, r_m = distancia a la Luna.
///
/// En el eje x (y = z = 0):
///
/// dΩ/dx = x - (1-μ)(x+μ)/r_e³ - μ(x-1+μ)/r_m³ = 0
///
/// Semilla inicial: aproximación asintótica de Richardson de 3er orden.
pub fn l1_position() -> f64 {
    let mu = MU;
    let xi0 = (mu / 3.0).cbrt() * (1.0 - mu / 3.0 - mu * mu / 9.0 - 23.0 * mu * mu * mu / 81.0);
    let mut x = 1.0 - mu - xi0;

    for _ in 0..20 {
        let r1 = x + mu;
        let r2 = 1.0 - mu - x;
        let r1_3 = r1 * r1 * r1;
        let r2_3 = r2 * r2 * r2;

        let f = x - (1.0 - mu) * r1 / r1_3 - mu * (x - 1.0 + mu) / r2_3;
        let df = 1.0 + 2.0 * (1.0 - mu) / r1_3 + 2.0 * mu / r2_3;

        let dx = f / df;
        x -= dx;

        if dx.abs() < 1e-15 {
            break;
        }
    }
    x
}

// ============================================================================
// DERIVADAS SEGUNDAS DEL PSEUDO-POTENCIAL EN L1
// ============================================================================

/// Calcula las derivadas segundas del pseudo-potencial Ω en L1 = (x_L1, 0, 0).
///
/// Retorna (Ω_xx, Ω_yy, Ω_zz).
///
/// Evaluado en y = 0, z = 0:
///   Ω_xx = 1 + 4(1-μ)/r_e³ + 4μ/r_m³
///   Ω_yy = 1 - (1-μ)/r_e³ - μ/r_m³
///   Ω_zz = -(1-μ)/r_e³ - μ/r_m³
fn omega_second_derivatives(xl1: f64) -> (f64, f64, f64) {
    let mu = MU;
    let r1 = xl1 + mu;
    let r2 = 1.0 - mu - xl1;
    let r1_3 = r1 * r1 * r1;
    let r2_3 = r2 * r2 * r2;

    let omega_xx = 1.0 + 4.0 * (1.0 - mu) / r1_3 + 4.0 * mu / r2_3;
    let omega_yy = 1.0 - (1.0 - mu) / r1_3 - mu / r2_3;
    let omega_zz = -(1.0 - mu) / r1_3 - mu / r2_3;

    (omega_xx, omega_yy, omega_zz)
}

// ============================================================================
// AUTOVALOR INESTABLE DE L1
// ============================================================================

/// Calcula el autovalor inestable (real positivo) de la linealización en L1.
///
/// De la ecuación característica del subsistema (x, y, vx, vy):
///   λ⁴ + (4 - a)λ² - (a² - 4) = 0
/// donde a = 1 + 2Ω_xx.
///
/// La raíz positiva es: λ_u = sqrt((a + sqrt(a² - 4)) / 2)
pub fn unstable_eigenvalue() -> f64 {
    let xl1 = l1_position();
    let (omega_xx, _, _) = omega_second_derivatives(xl1);
    let a = 1.0 + 2.0 * omega_xx;
    let discriminant = a * a - 4.0;
    let lambda_sq = (a + discriminant.sqrt()) / 2.0;
    lambda_sq.sqrt()
}

// ============================================================================
// FACTOR DE FORMA DEL AUTOVECTOR
// ============================================================================

/// Calcula el factor de forma k = δy/δx para los autovectores en L1.
///
/// De la ecuación ÿ - 2ẋ = Ω_yy · y, con la forma del autovector:
///   k(λ² - Ω_yy) = 2λ  →  k = 2λ / (λ² - Ω_yy)
///
/// Para el sistema Tierra-Luna en L1: k ≈ 0.476 (positivo).
pub fn eigenvector_y_factor() -> f64 {
    let lambda = unstable_eigenvalue();
    let xl1 = l1_position();
    let (_, omega_yy, _) = omega_second_derivatives(xl1);
    let denominator = lambda * lambda - omega_yy;
    2.0 * lambda / denominator
}

// ============================================================================
// AUTOVECTORES COMPLETOS
// ============================================================================

/// Construye el autovector inestable normalizado con δx = 1.
/// v_inestable = (1, k, 0, λ, λk, 0)
pub fn unstable_eigenvector() -> [f64; 6] {
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();
    [1.0, k, 0.0, lambda, lambda * k, 0.0]
}

/// Construye el autovector estable normalizado con δx = 1.
/// v_estable = (1, k, 0, -λ, -λk, 0)
pub fn stable_eigenvector() -> [f64; 6] {
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();
    [1.0, k, 0.0, -lambda, -lambda * k, 0.0]
}

// ============================================================================
// SENSIBILIDAD DE LA CONSTANTE DE JACOBI
// ============================================================================

/// Constante geométrica del sistema: sensibilidad de C_J respecto a ε² en L1.
///
/// Para cualquier dirección basada en los autovectores normalizados
/// (δx = ±1, δy = ±k, δvx = ±λ, δvy = ±λk), el valor de α es idéntico
/// porque los signos se cancelan al elevar al cuadrado:
///
///   α = Ω_xx + k²·Ω_yy - λ²(1 + k²)
///
/// Esta constante depende solo de μ (el sistema planetario), no de la
/// dirección elegida ni de target_cj. Se precalcula una vez y se reutiliza
/// en todos los ajustes de Newton-Raphson sobre ε².
///
/// # Uso
/// ```
/// let alpha = jacobi_sensitivity_l1();
/// let eps_sq_new = eps_sq - (cj_current - target_cj) / alpha;
/// ```
pub fn jacobi_sensitivity_l1() -> f64 {
    let xl1 = l1_position();
    let (omega_xx, omega_yy, _) = omega_second_derivatives(xl1);
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();
    omega_xx + k * k * omega_yy - lambda * lambda * (1.0 + k * k)
}

/// Calcula la derivada segunda de C_J a lo largo de una dirección arbitraria.
///
/// Forma cuadrática exacta: α = (Ω_xx·δx² + Ω_yy·δy² + Ω_zz·δz²) 
///                               - (δv_x² + δv_y² + δv_z²)
///
/// Para direcciones basadas en autovectores, `jacobi_sensitivity_l1()` es
/// equivalente y más eficiente (no requiere argumento).
pub fn jacobi_sensitivity_along_eigenvector(direction: &[f64; 6]) -> f64 {
    let xl1 = l1_position();
    let (omega_xx, omega_yy, omega_zz) = omega_second_derivatives(xl1);

    let spatial_part = omega_xx * direction[0] * direction[0]
                     + omega_yy * direction[1] * direction[1]
                     + omega_zz * direction[2] * direction[2];
    let velocity_part = direction[3] * direction[3]
                      + direction[4] * direction[4]
                      + direction[5] * direction[5];

    spatial_part - velocity_part
}

// ============================================================================
// TIPOS DE ERROR
// ============================================================================

#[derive(Error, Debug)]
pub enum ManifoldError {
    #[error("Error en conexión: {0}")]
    ConnectionError(String),
    #[error("No se pudo alcanzar la constante de Jacobi objetivo: C_J={0:.6}, objetivo={1:.6}")]
    JacobiConvergenceError(f64, f64),
}

// ============================================================================
// TIPOS DE VARIEDAD
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ManifoldType {
    /// Variedad estable: converge a L1 en tiempo positivo (t → +∞).
    Stable,
    /// Variedad inestable: diverge de L1 en tiempo positivo (t → +∞).
    Unstable,
    /// Tránsito Tierra → Luna: variedad estable propagada hacia adelante.
    /// La trayectoria converge a L1 desde la cuenca terrestre, lo cruza,
    /// y emerge en la cuenca lunar por la variedad inestable.
    TransitToMoon,
    /// Tránsito Luna → Tierra: variedad inestable propagada hacia atrás.
    /// La trayectoria converge a L1 desde la cuenca lunar, lo cruza,
    /// y emerge en la cuenca terrestre por la variedad estable.
    TransitToEarth,
}

// ============================================================================
// GENERACIÓN DE PUNTOS EN LA VARIEDAD ESTABLE (RAMA TERRESTRE)
// ============================================================================

/// Genera un punto sobre la variedad ESTABLE de L1, rama terrestre.
///
/// **Rama terrestre (inyección Tierra → L1):**
///   - δx < 0 (hacia la Tierra, x < x_L1)
///   - δy < 0 (y negativo)
///   - vx > 0 (velocidad hacia L1 en tiempo positivo)
///   - vy > 0
///   - Dirección: (-1, -k, 0, λ, λk, 0)
///
/// **Ajuste de ε:** Newton sobre ε² usando α = `jacobi_sensitivity_l1()`.
pub fn stable_manifold_point_terrestrial(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    let xl1 = l1_position();
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();

    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    if target_cj <= cj_l1 + 1e-12 {
        return Err(ManifoldError::JacobiConvergenceError(target_cj, cj_l1));
    }

    // Rama terrestre: δx < 0, δy < 0, vx > 0 (hacia L1)
    let direction = [-1.0, -k, 0.0, lambda, lambda * k, 0.0];
    let alpha = jacobi_sensitivity_l1();

    let mut eps = amplitude_guess;
    for _ in 0..30 {
        let state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj = jacobi_constant(&state);
        let err = cj - target_cj;

        if err.abs() < 1e-10 {
            return Ok(state);
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            return Err(ManifoldError::JacobiConvergenceError(cj, target_cj));
        }

        eps = eps_sq_new.sqrt();
        eps = eps.clamp(0.0001, 0.2);
    }

    Err(ManifoldError::JacobiConvergenceError(
        jacobi_constant(&StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            0.0,
            eps * direction[3],
            eps * direction[4],
            0.0,
        )),
        target_cj,
    ))
}

// ============================================================================
// GENERACIÓN DE PUNTOS EN LA VARIEDAD ESTABLE (RAMA LUNAR)
// ============================================================================

/// Genera un punto sobre la variedad ESTABLE de L1, rama lunar.
///
/// **Rama lunar (P₃ en el artículo):**
///   - δx < 0 (hacia la Luna desde L1, x < x_L1)
///   - δy < 0 (y negativo)
///   - vx > 0, vy > 0
///   - Dirección: (-1, -k, 0, λ, λk, 0)
///   - Integrando hacia atrás en el tiempo se obtiene la variedad que llega
///     a la Luna desde L1.
///
/// Nota: Esta función es idéntica a `stable_manifold_point_terrestrial` en
/// su implementación actual. Se mantiene como alias semántico. La diferencia
/// entre rama terrestre y lunar está en cómo se propaga el estado (hacia
/// adelante o hacia atrás en el tiempo), no en la construcción del punto.
pub fn stable_manifold_point_lunar(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    stable_manifold_point_terrestrial(target_cj, amplitude_guess)
}

/// Alias mantenido por compatibilidad con código existente.
#[deprecated(since = "2.4.0", note = "Usar `stable_manifold_point_terrestrial` o `stable_manifold_point_lunar`")]
pub fn stable_manifold_point_for_jacobi(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    stable_manifold_point_terrestrial(target_cj, amplitude_guess)
}

// ============================================================================
// PUNTO DE TRÁNSITO LUNA → TIERRA
// ============================================================================

/// Genera un punto para tránsito Luna → Tierra a través de L1.
///
/// **Rama de retorno (desde L1 hacia la Tierra):**
///   - δx < 0 (lado terrestre, x < x_L1)
///   - δy < 0
///   - vx < 0, vy < 0 (alejándose de L1 hacia Tierra)
///   - Dirección: (-1, -k, 0, -λ, -λk, 0)
///
/// Se propaga hacia atrás en el tiempo (time_direction = -1).
pub fn transit_point_to_earth(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    let xl1 = l1_position();
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();

    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    if target_cj >= cj_l1 + 1e-12 {
        return Err(ManifoldError::JacobiConvergenceError(target_cj, cj_l1));
    }

    // Lado terrestre, alejándose hacia Tierra: δx < 0, δy < 0, vx < 0, vy < 0
    let direction = [-1.0, -k, 0.0, -lambda, -lambda * k, 0.0];
    let alpha = jacobi_sensitivity_l1();

    let mut eps = amplitude_guess;
    for _ in 0..30 {
        let state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj = jacobi_constant(&state);
        let err = cj - target_cj;

        if err.abs() < 1e-10 {
            return Ok(state);
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            return Err(ManifoldError::JacobiConvergenceError(cj, target_cj));
        }

        eps = eps_sq_new.sqrt();
        eps = eps.clamp(0.0001, 0.2);
    }

    Err(ManifoldError::JacobiConvergenceError(
        jacobi_constant(&StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            0.0,
            eps * direction[3],
            eps * direction[4],
            0.0,
        )),
        target_cj,
    ))
}

// ============================================================================
// PUNTO DE TRÁNSITO TIERRA → LUNA
// ============================================================================

/// Genera un punto para tránsito Tierra → Luna a través de L1.
///
/// **Rama de tránsito (desde L1 hacia la Luna):**
///   - δx > 0 (lado lunar, x > x_L1)
///   - δy > 0
///   - vx > 0, vy > 0 (alejándose de L1)
///   - Dirección: (1, k, 0, λ, λk, 0)
///
/// Se propaga hacia adelante en el tiempo (time_direction = +1).
pub fn transit_point_to_moon(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    let xl1 = l1_position();
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();

    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    if target_cj >= cj_l1 + 1e-12 {
        return Err(ManifoldError::JacobiConvergenceError(target_cj, cj_l1));
    }

    // Lado lunar, alejándose: δx > 0, δy > 0, vx > 0, vy > 0
    let direction = [1.0, k, 0.0, lambda, lambda * k, 0.0];
    let alpha = jacobi_sensitivity_l1();

    let mut eps = amplitude_guess;
    for _ in 0..30 {
        let state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj = jacobi_constant(&state);
        let err = cj - target_cj;

        if err.abs() < 1e-10 {
            return Ok(state);
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            return Err(ManifoldError::JacobiConvergenceError(cj, target_cj));
        }

        eps = eps_sq_new.sqrt();
        eps = eps.clamp(0.0001, 0.2);
    }

    Err(ManifoldError::JacobiConvergenceError(
        jacobi_constant(&StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            0.0,
            eps * direction[3],
            eps * direction[4],
            0.0,
        )),
        target_cj,
    ))
}


// ============================================================================
// GENERACIÓN DE PUNTOS EN LA VARIEDAD INESTABLE
// ============================================================================

/// Genera un punto sobre la variedad INESTABLE de L1, rama exterior.
///
/// **Rama exterior (escape de L1 hacia la Tierra):**
///   - δx > 0, δy > 0, vx > 0, vy > 0
///   - Dirección: (1, k, 0, λ, λk, 0)
///   - Integrando hacia adelante, la trayectoria se aleja de L1 hacia la Tierra.
pub fn unstable_manifold_point_exterior(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    let xl1 = l1_position();
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();

    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    if target_cj <= cj_l1 + 1e-12 {
        return Err(ManifoldError::JacobiConvergenceError(target_cj, cj_l1));
    }

    let direction = [1.0, k, 0.0, lambda, -lambda * k, 0.0];
    let alpha = jacobi_sensitivity_l1();

    let mut eps = amplitude_guess;
    for _ in 0..30 {
        let state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj = jacobi_constant(&state);
        let err = cj - target_cj;

        if err.abs() < 1e-10 {
            return Ok(state);
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            return Err(ManifoldError::JacobiConvergenceError(cj, target_cj));
        }

        eps = eps_sq_new.sqrt();
        eps = eps.clamp(0.0001, 0.2);
    }

    Err(ManifoldError::JacobiConvergenceError(
        jacobi_constant(&StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            0.0,
            eps * direction[3],
            eps * direction[4],
            0.0,
        )),
        target_cj,
    ))
}

/// Alias mantenido por compatibilidad.
#[deprecated(since = "2.4.0", note = "Usar `unstable_manifold_point_exterior`")]
pub fn unstable_manifold_point_for_jacobi(
    target_cj: f64,
    amplitude_guess: f64,
) -> Result<StateVector, ManifoldError> {
    unstable_manifold_point_exterior(target_cj, amplitude_guess)
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_l1_position_consistency() {
        let xl1 = l1_position();
        let mu = MU;
        assert!(xl1 > -mu && xl1 < 1.0 - mu);

        let r1 = xl1 + mu;
        let r2 = 1.0 - mu - xl1;
        let grad = xl1 - (1.0 - mu) * r1 / (r1 * r1 * r1) - mu * (xl1 - 1.0 + mu) / (r2 * r2 * r2);
        assert!(grad.abs() < 1e-12, "Gradiente en L1 no es cero: {}", grad);
    }

    #[test]
    fn test_eigenvector_factor_derivation() {
        let lambda = unstable_eigenvalue();
        let xl1 = l1_position();
        let (_, omega_yy, _) = omega_second_derivatives(xl1);
        let k = eigenvector_y_factor();
        let lhs = lambda * lambda * k - 2.0 * lambda;
        let rhs = omega_yy * k;
        assert!((lhs - rhs).abs() < 1e-12,
            "Factor k no satisface ecuación: {} ≠ {}", lhs, rhs);
    }

    #[test]
    fn test_sensitivity_l1_matches_directional() {
        let lambda = unstable_eigenvalue();
        let k = eigenvector_y_factor();
        let alpha_direct = jacobi_sensitivity_l1();

        // Debe coincidir con cualquier dirección de autovector
        let dir1 = [-1.0, -k, 0.0, lambda, lambda * k, 0.0];
        let dir2 = [1.0, k, 0.0, lambda, lambda * k, 0.0];
        let dir3 = [1.0, k, 0.0, -lambda, -lambda * k, 0.0];

        let alpha1 = jacobi_sensitivity_along_eigenvector(&dir1);
        let alpha2 = jacobi_sensitivity_along_eigenvector(&dir2);
        let alpha3 = jacobi_sensitivity_along_eigenvector(&dir3);

        assert!((alpha_direct - alpha1).abs() < 1e-15);
        assert!((alpha_direct - alpha2).abs() < 1e-15);
        assert!((alpha_direct - alpha3).abs() < 1e-15);
    }

    #[test]
    fn test_analytical_vs_finite_difference_sensitivity() {
        let xl1 = l1_position();
        let lambda = unstable_eigenvalue();
        let k = eigenvector_y_factor();
        let direction = [-1.0, -k, 0.0, lambda, lambda * k, 0.0];

        let alpha_analytic = jacobi_sensitivity_l1();

        let eps = 0.001;
        let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj_l1 = jacobi_constant(&state_l1);
        let state_pert = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            0.0,
            eps * direction[3],
            eps * direction[4],
            0.0,
        );
        let alpha_fd = (jacobi_constant(&state_pert) - cj_l1) / (eps * eps);

        let relative_error = (alpha_analytic - alpha_fd).abs() / alpha_analytic.abs();
        assert!(relative_error < 1e-6,
            "α analítica={} vs α fd={}, error_rel={}", alpha_analytic, alpha_fd, relative_error);
    }

    #[test]
    fn test_stable_manifold_terrestrial() {
        let xl1 = l1_position();
        let cj_l1 = jacobi_constant(&StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0));
        let target = cj_l1 + 0.01;

        let result = stable_manifold_point_terrestrial(target, 0.01);
        assert!(result.is_ok(), "Fallo: {:?}", result.err());

        let state = result.unwrap();
        let cj = jacobi_constant(&state);
        assert!((cj - target).abs() < 1e-8);
        assert!(state.x < xl1, "x={} debe ser < x_L1={}", state.x, xl1);
        assert!(state.y < 0.0, "y={} debe ser negativo", state.y);
        assert!(state.vx > 0.0, "vx={} debe ser positivo (hacia L1)", state.vx);
    }

    #[test]
    fn test_unstable_manifold_exterior() {
        let xl1 = l1_position();
        let cj_l1 = jacobi_constant(&StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0));
        let target = cj_l1 + 0.01;

        let result = unstable_manifold_point_exterior(target, 0.01);
        assert!(result.is_ok(), "Fallo: {:?}", result.err());

        let state = result.unwrap();
        let cj = jacobi_constant(&state);
        assert!((cj - target).abs() < 1e-8);
        assert!(state.x > xl1, "x={} debe ser > x_L1={}", state.x, xl1);
        assert!(state.y > 0.0, "y={} debe ser positivo", state.y);
    }

    #[test]
    fn test_jacobi_increases_with_amplitude() {
        let xl1 = l1_position();
        let cj_l1 = jacobi_constant(&StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0));
        let r1 = stable_manifold_point_terrestrial(cj_l1 + 0.005, 0.005).unwrap();
        let r2 = stable_manifold_point_terrestrial(cj_l1 + 0.010, 0.010).unwrap();
        assert!(jacobi_constant(&r1) < jacobi_constant(&r2));
    }

    #[test]
    fn test_rejects_invalid_jacobi() {
        let xl1 = l1_position();
        let cj_l1 = jacobi_constant(&StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert!(stable_manifold_point_terrestrial(cj_l1 - 0.01, 0.01).is_err());
        assert!(stable_manifold_point_terrestrial(cj_l1, 0.01).is_err());
    }

    #[test]
    fn test_eigenvector_properties() {
        let vu = unstable_eigenvector();
        let vs = stable_eigenvector();
        assert!((vu[1] - vs[1]).abs() < 1e-15, "δy deben ser iguales");
        assert!((vu[3] + vs[3]).abs() < 1e-15, "vx deben ser opuestos");
        assert!((vu[4] + vs[4]).abs() < 1e-15, "vy deben ser opuestos");
        assert!(vu[2].abs() < 1e-15 && vu[5].abs() < 1e-15);
        assert!(vs[2].abs() < 1e-15 && vs[5].abs() < 1e-15);
    }

    #[test]
    fn test_sensitivity_is_positive() {
        let alpha = jacobi_sensitivity_l1();
        assert!(alpha > 0.0, "α debe ser positivo: {}", alpha);
    }
}