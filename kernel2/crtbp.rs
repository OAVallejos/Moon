//! crtbp.rs Dinámica del Problema Restringido de Tres Cuerpos Circular (CRTBP)                    
//! Ecuaciones 5-6 del artículo Almeida Jr. et al. (2026)   
//!
//! Sistema de ecuaciones adimensional en marco rotante sinódico
//!
//! UNIDADES:
//! - Longitud: Unidades Canónicas (LU) donde 1 LU = D_CHAR = 384,400 km
//! - Tiempo: Unidades Canónicas (TU) donde 1 TU = 1/ω ≈ 4.342 días
//! - Velocidad: Unidades Canónicas (VU) donde 1 VU = D_CHAR * ω ≈ 1.0245 km/s
//!
//! La constante de Jacobi C_J debe ser ~3.20034 en L1 para unidades canónicas.

use nalgebra::SVector;

/// Vector de estado 6D: [x, y, z, vx, vy, vz] (adimensional)
pub type StateVector = SVector<f64, 6>;
/// Vector de derivadas 6D: [vx, vy, vz, ax, ay, az]
pub type DerivativeVector = SVector<f64, 6>;

/// Parámetro de masa normalizado μ = M_Luna / (M_Tierra + M_Luna)
pub const MU: f64 = 0.01215058560962404;

/// Posición de la Tierra en el eje x (primario mayor)
pub const X_EARTH: f64 = -MU;
/// Posición de la Luna en el eje x (primario menor)
pub const X_MOON: f64 = 1.0 - MU;

/// Distancia al punto L1 (colineal, entre Tierra y Luna)
pub const X_L1: f64 = 0.8369151948720569;

/// Constante de Jacobi teórica en L1 (unidades canónicas)
pub const CJ_L1_THEORETICAL: f64 = 3.200344;

// ============================================================================
// CONSTANTES DE UNIDADES PARA CONVERSIÓN
// ============================================================================

/// Distancia característica [km] (Tierra-Luna)
pub const D_CHAR_KM: f64 = 384_400.0;

/// Velocidad característica [km/s]
pub const V_CHAR_KMS: f64 = 1.02458;

// ============================================================================
// DERIVADAS DEL CRTBP
// ============================================================================

/// Evalúa las derivadas del CRTBP (Ecuación 5 del artículo)
///
/// # Argumentos
/// * `_t` - Tiempo adimensional (no autónomo, no se usa)
/// * `state` - Vector de estado [x, y, z, vx, vy, vz] en unidades canónicas
///
/// # Devuelve
/// Vector de derivadas [vx, vy, vz, ax, ay, az]
pub fn crtbp_derivatives(_t: f64, state: &StateVector) -> DerivativeVector {
    let x = state[0];
    let y = state[1];
    let z = state[2];
    let vx = state[3];
    let vy = state[4];
    let vz = state[5];

    // Distancias a los primarios
    let r_e = distance_to_earth(x, y, z);
    let r_m = distance_to_moon(x, y, z);

    // Términos gravitacionales
    let factor_e = (1.0 - MU) / (r_e * r_e * r_e);
    let factor_m = MU / (r_m * r_m * r_m);

    // Pseudo-potencial efectivo y aceleraciones
    // Ω = (1/2)(x² + y²) + (1-μ)/r_e + μ/r_m
    // Ecuaciones de movimiento en marco rotante:
    // ẍ - 2ẏ = ∂Ω/∂x
    // ÿ + 2ẋ = ∂Ω/∂y
    // z̈ = ∂Ω/∂z

    let domega_dx = x - factor_e * (x - X_EARTH) - factor_m * (x - X_MOON);
    let domega_dy = y - factor_e * y - factor_m * y;
    let domega_dz = -factor_e * z - factor_m * z;

    // Términos de Coriolis
    let ax = 2.0 * vy + domega_dx;
    let ay = -2.0 * vx + domega_dy;
    let az = domega_dz;

    DerivativeVector::new(vx, vy, vz, ax, ay, az)
}

// ============================================================================
// CONSTANTE DE JACOBI (CORREGIDA CON UNIDADES)
// ============================================================================

/// Calcula la constante de Jacobi (Ecuación 6 del artículo)
/// C_J = 2Ω - (vx² + vy² + vz²)
/// donde Ω es el pseudo-potencial efectivo
///
/// # IMPORTANTE: Unidades
/// El estado DEBE estar en unidades canónicas (LU, VU).
/// Si el estado está en km y km/s, se convierte automáticamente.
///
/// # Valores esperados
/// - En L1: C_J ≈ 3.200344
/// - Para cruzar L1: C_J < 3.200344
/// - Para variedades estables/inestables: C_J > 3.200344
pub fn jacobi_constant(state: &StateVector) -> f64 {
    // Detectar unidades y convertir si es necesario
    // Si |x| > 10, asumimos que está en km (ya que 1 LU ≈ 384,400 km)
    let (x, y, z, vx, vy, vz) = if state[0].abs() > 10.0 {
        // Convertir de km y km/s a unidades canónicas
        (
            state[0] / D_CHAR_KM,
            state[1] / D_CHAR_KM,
            state[2] / D_CHAR_KM,
            state[3] / V_CHAR_KMS,
            state[4] / V_CHAR_KMS,
            state[5] / V_CHAR_KMS,
        )
    } else {
        // Ya está en unidades canónicas
        (state[0], state[1], state[2], state[3], state[4], state[5])
    };

    let r_e = distance_to_earth(x, y, z);
    let r_m = distance_to_moon(x, y, z);

    // Pseudo-potencial: Ω = (1/2)(x² + y²) + (1-μ)/r_e + μ/r_m
    let omega = 0.5 * (x * x + y * y) + (1.0 - MU) / r_e + MU / r_m;

    // C_J = 2Ω - v²
    2.0 * omega - (vx * vx + vy * vy + vz * vz)
}

/// Versión de jacobi_constant que acepta un array de 6 elementos
pub fn jacobi_constant_array(state: &[f64; 6]) -> f64 {
    let sv = StateVector::new(state[0], state[1], state[2], state[3], state[4], state[5]);
    jacobi_constant(&sv)
}

/// Versión de jacobi_constant que acepta un slice
pub fn jacobi_constant_slice(state: &[f64]) -> f64 {
    if state.len() < 6 {
        return f64::NAN;
    }
    let sv = StateVector::new(state[0], state[1], state[2], state[3], state[4], state[5]);
    jacobi_constant(&sv)
}

/// Verifica si un estado puede cruzar L1 (C_J < C_J(L1))
pub fn can_cross_l1(state: &StateVector) -> bool {
    jacobi_constant(state) < CJ_L1_THEORETICAL
}

/// Verifica si un estado está en una variedad (C_J > C_J(L1))
pub fn is_on_manifold(state: &StateVector) -> bool {
    jacobi_constant(state) > CJ_L1_THEORETICAL
}

// ============================================================================
// FUNCIONES DE DISTANCIA
// ============================================================================

/// Distancia a la Tierra (primario mayor)
pub fn distance_to_earth(x: f64, y: f64, z: f64) -> f64 {
    let dx = x - X_EARTH;
    (dx * dx + y * y + z * z).sqrt()
}

/// Distancia a la Luna (primario menor)
pub fn distance_to_moon(x: f64, y: f64, z: f64) -> f64 {
    let dx = x - X_MOON;
    (dx * dx + y * y + z * z).sqrt()
}

// ============================================================================
// CONVERSIONES DE COORDENADAS
// ============================================================================

/// Convierte de coordenadas polares geocéntricas (r, θ) a cartesianas sinódicas
/// Ecuaciones 13-14 del artículo
pub fn polar_to_cartesian(r: f64, theta: f64, dr: f64, dtheta: f64) -> StateVector {
    let x = r * theta.cos() - MU;
    let y = r * theta.sin();
    let z = 0.0;
    let vx = dr * theta.cos() - r * dtheta * theta.sin();
    let vy = dr * theta.sin() + r * dtheta * theta.cos();
    let vz = 0.0;
    StateVector::new(x, y, z, vx, vy, vz)
}

/// Convierte de coordenadas polares selenocéntricas (R, Θ) a cartesianas sinódicas
/// Ecuaciones 21-22 del artículo
pub fn lunar_polar_to_cartesian(r: f64, theta: f64, dr: f64, dtheta: f64) -> StateVector {
    let x = r * theta.cos() + (1.0 - MU);
    let y = r * theta.sin();
    let z = 0.0;
    let vx = dr * theta.cos() - r * dtheta * theta.sin();
    let vy = dr * theta.sin() + r * dtheta * theta.cos();
    let vz = 0.0;
    StateVector::new(x, y, z, vx, vy, vz)
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jacobi_at_l1() {
        // En el punto L1, la velocidad es cero
        let state_l1 = StateVector::new(X_L1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj = jacobi_constant(&state_l1);
        // La constante de Jacobi en L1 debe ser ~3.20034 (valor conocido)
        let diff = (cj - CJ_L1_THEORETICAL).abs();
        assert!(diff < 0.001, "C_J(L1) = {:.8}, esperado {:.8}", cj, CJ_L1_THEORETICAL);
        println!("✅ C_J(L1) = {:.8} (teórico: {:.8})", cj, CJ_L1_THEORETICAL);
    }

    #[test]
    fn test_jacobi_with_km_units() {
        // Estado en km y km/s (debe convertirse automáticamente)
        let xl1_km = X_L1 * D_CHAR_KM;
        let state_km = StateVector::new(xl1_km, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj = jacobi_constant(&state_km);
        let diff = (cj - CJ_L1_THEORETICAL).abs();
        assert!(diff < 0.001, "C_J desde km = {:.8}, esperado {:.8}", cj, CJ_L1_THEORETICAL);
        println!("✅ C_J desde km = {:.8} (teórico: {:.8})", cj, CJ_L1_THEORETICAL);
    }

    #[test]
    fn test_can_cross_l1() {
        // Estado con C_J < C_J(L1) - puede cruzar
        let state = StateVector::new(X_L1 - 0.01, -0.001, 0.0, 0.01, 0.005, 0.0);
        let cj = jacobi_constant(&state);
        assert!(cj < CJ_L1_THEORETICAL, "C_J = {:.8} debe ser < {:.8}", cj, CJ_L1_THEORETICAL);
        assert!(can_cross_l1(&state));
        println!("✅ C_J = {:.8} < {:.8} → puede cruzar L1", cj, CJ_L1_THEORETICAL);
    }

    #[test]
    fn test_is_on_manifold() {
        // Estado con C_J > C_J(L1) - está en variedad
        let state = StateVector::new(X_L1 + 0.01, 0.001, 0.0, 0.01, 0.005, 0.0);
        let cj = jacobi_constant(&state);
        assert!(cj > CJ_L1_THEORETICAL, "C_J = {:.8} debe ser > {:.8}", cj, CJ_L1_THEORETICAL);
        assert!(is_on_manifold(&state));
        println!("✅ C_J = {:.8} > {:.8} → está en variedad", cj, CJ_L1_THEORETICAL);
    }

    #[test]
    fn test_derivatives_at_l1_equilibrium() {
        // En L1 con velocidad cero, las derivadas deben ser cero (equilibrio)
        let state_l1 = StateVector::new(X_L1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let deriv = crtbp_derivatives(0.0, &state_l1);
        // La aceleración debe ser aproximadamente cero en el punto de equilibrio
        assert!(deriv[3].abs() < 1e-10);
        assert!(deriv[4].abs() < 1e-10);
        println!("✅ Derivadas en L1: ax={:.2e}, ay={:.2e}", deriv[3], deriv[4]);
    }

    #[test]
    fn test_jacobi_array_version() {
        let state = [X_L1, 0.0, 0.0, 0.0, 0.0, 0.0];
        let cj = jacobi_constant_array(&state);
        let diff = (cj - CJ_L1_THEORETICAL).abs();
        assert!(diff < 0.001);
        println!("✅ jacobi_constant_array = {:.8}", cj);
    }

    #[test]
    fn test_jacobi_slice_version() {
        let state = vec![X_L1, 0.0, 0.0, 0.0, 0.0, 0.0];
        let cj = jacobi_constant_slice(&state);
        let diff = (cj - CJ_L1_THEORETICAL).abs();
        assert!(diff < 0.001);
        println!("✅ jacobi_constant_slice = {:.8}", cj);
    }
}
