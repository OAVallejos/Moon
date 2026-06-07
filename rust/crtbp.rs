//! Dinámica del Problema Restringido de Tres Cuerpos Circular (CRTBP)
//! Ecuaciones 5-6 del artículo Almeida Jr. et al. (2026)
//! 
//! Sistema de ecuaciones adimensional en marco rotante sinódico

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

/// Evalúa las derivadas del CRTBP (Ecuación 5 del artículo)
/// 
/// # Argumentos
/// * `_t` - Tiempo adimensional (no autónomo, no se usa)
/// * `state` - Vector de estado [x, y, z, vx, vy, vz]
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

/// Calcula la constante de Jacobi (Ecuación 6 del artículo)
/// C_J = 2Ω - (vx² + vy² + vz²)
/// donde Ω es el pseudo-potencial efectivo
pub fn jacobi_constant(state: &StateVector) -> f64 {
    let x = state[0];
    let y = state[1];
    let z = state[2];
    let vx = state[3];
    let vy = state[4];
    let vz = state[5];
    
    let r_e = distance_to_earth(x, y, z);
    let r_m = distance_to_moon(x, y, z);
    
    // Pseudo-potencial: Ω = (1/2)(x² + y²) + (1-μ)/r_e + μ/r_m
    let omega = 0.5 * (x * x + y * y) + (1.0 - MU) / r_e + MU / r_m;
    
    // C_J = 2Ω - v²
    2.0 * omega - (vx * vx + vy * vy + vz * vz)
}

/// Distancia a la Tierra (primario mayor)
fn distance_to_earth(x: f64, y: f64, z: f64) -> f64 {
    let dx = x - X_EARTH;
    (dx * dx + y * y + z * z).sqrt()
}

/// Distancia a la Luna (primario menor)
fn distance_to_moon(x: f64, y: f64, z: f64) -> f64 {
    let dx = x - X_MOON;
    (dx * dx + y * y + z * z).sqrt()
}

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

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_jacobi_at_l1() {
        // En el punto L1, la velocidad es cero
        let state_l1 = StateVector::new(X_L1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj = jacobi_constant(&state_l1);
        // La constante de Jacobi en L1 debe ser ~3.20034 (valor conocido)
        assert!((cj - 3.20034).abs() < 0.001);
    }
    
    #[test]
    fn test_derivatives_at_l1_equilibrium() {
        // En L1 con velocidad cero, las derivadas deben ser cero (equilibrio)
        let state_l1 = StateVector::new(X_L1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let deriv = crtbp_derivatives(0.0, &state_l1);
        // La aceleración debe ser aproximadamente cero en el punto de equilibrio
        assert!(deriv[3].abs() < 1e-10);
        assert!(deriv[4].abs() < 1e-10);
    }
}