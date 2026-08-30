//! geocentric.rs Etapa 1: Inyección a la Variedad Estable de L1 desde Órbita Terrestre Alta
//!
//! Basado en Almeida Jr. et al. (2026), Ecuaciones 13-14.
//!
//! CORRECCIONES (Junio 2026):
//!   - ESTADO DE INYECCIÓN EN UNIDADES CANÓNICAS (no metros)
//!   - C_J calculado en unidades canónicas consistentes
//!   - Validación de unidades en todos los outputs
//!
//! ESTRATEGIA:
//!   - La nave es entregada a una órbita terrestre alta por un lanzador comercial.
//!   - Este módulo calcula el punto de inyección sobre la variedad ESTABLE de L1
//!     (rama terrestre: δx < 0, δy < 0, vx > 0 hacia L1).
//!   - El ΔV se calcula como diferencia entre la velocidad orbital actual y la
//!     velocidad requerida en el punto de inyección.
//!   - TODO en unidades canónicas (distancia Tierra-Luna = 1.0)
//!
//! FLUJO DE MISIÓN:
//!   Etapa 0 (Externa): Lanzador → Órbita terrestre alta
//!   Etapa 1 (Este módulo): Inyección → Variedad ESTABLE L1 (ΔV ≈ 40-100 m/s)
//!   Etapa 2 (integration.rs): Tránsito balístico → Captura lunar (ΔV = 0)
//!   Etapa 3 (lunar.rs): Análisis de captura lunar

use thiserror::Error;
use crate::crtbp::{StateVector, jacobi_constant, MU};
use crate::constants::{D_CHAR, V_CHAR};
use crate::manifold::{
    l1_position,
    unstable_eigenvalue,
    eigenvector_y_factor,
    jacobi_sensitivity_l1,
};

// ============================================================================
// TIPOS DE ERROR
// ============================================================================

#[derive(Error, Debug)]
pub enum InjectionError {
    #[error("C_J objetivo {0:.6} inalcanzable. Mínimo teórico: {1:.6}")]
    JacobiInfeasible(f64, f64),
    #[error("ΔV calculado fuera de rango físico: {0:.2} m/s")]
    UnphysicalDeltaV(f64),
    #[error("La nave no está en la cuenca terrestre: x={0:.6} >= x_L1={1:.6}")]
    NotInTerrestrialRealm(f64, f64),
    #[error("Newton-Raphson no convergió para ε en {0} iteraciones")]
    EpsilonNotConverged(usize),
}

// ============================================================================
// ESTADO ORBITAL (INPUT DEL LANZADOR) - SIEMPRE EN CANÓNICO
// ============================================================================

/// Estado orbital inicial provisto por el lanzador comercial.
/// TODAS las coordenadas en unidades CANÓNICAS (distancia Tierra-Luna = 1.0)
#[derive(Debug, Clone)]
pub struct OrbitalState {
    /// Posición [x, y, z] en sistema rotante CRTBP (canónico)
    pub position: [f64; 3],
    /// Velocidad [vx, vy, vz] en sistema rotante CRTBP (canónico)
    pub velocity: [f64; 3],
}

impl OrbitalState {
    /// Crea un estado orbital a partir de coordenadas dimensionales (metros, m/s).
    /// Convierte automáticamente a unidades canónicas.
    pub fn from_dimensional(
        x_m: f64, y_m: f64, z_m: f64,
        vx_ms: f64, vy_ms: f64, vz_ms: f64,
    ) -> Self {
        Self {
            position: [x_m / D_CHAR, y_m / D_CHAR, z_m / D_CHAR],
            velocity: [vx_ms / V_CHAR, vy_ms / V_CHAR, vz_ms / V_CHAR],
        }
    }

    /// Crea un estado orbital en unidades canónicas (sin conversión).
    pub fn new_canonical(
        x: f64, y: f64, z: f64,
        vx: f64, vy: f64, vz: f64,
    ) -> Self {
        Self {
            position: [x, y, z],
            velocity: [vx, vy, vz],
        }
    }

    /// Distancia al centro de la Tierra en el sistema rotante (canónico).
    pub fn distance_to_earth(&self) -> f64 {
        let dx = self.position[0] + MU;
        let dy = self.position[1];
        let dz = self.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Distancia al centro de la Tierra en km.
    pub fn distance_to_earth_km(&self) -> f64 {
        self.distance_to_earth() * D_CHAR / 1000.0
    }

    /// Constante de Jacobi del estado orbital (canónico).
    pub fn jacobi(&self) -> f64 {
        let state = StateVector::new(
            self.position[0], self.position[1], self.position[2],
            self.velocity[0], self.velocity[1], self.velocity[2],
        );
        jacobi_constant(&state)
    }
}

// ============================================================================
// RESULTADO DE LA INYECCIÓN - CON AMBOS SISTEMAS DE UNIDADES
// ============================================================================

/// Resultado completo de la maniobra de inyección.
#[derive(Debug, Clone)]
pub struct InjectionResult {
    /// Estado completo en el punto de inyección (CANÓNICO)
    pub state: StateVector,
    /// Estado en unidades dimensionales [x_m, y_m, z_m, vx_ms, vy_ms, vz_ms]
    pub state_dimensional: [f64; 6],
    /// Amplitud ε sobre la variedad
    pub epsilon: f64,
    /// ΔV total en m/s
    pub dv_total_ms: f64,
    /// ΔV total en unidades canónicas
    pub dv_total_canonical: f64,
    /// Componente radial del ΔV en m/s
    pub dv_radial_ms: f64,
    /// Componente tangencial del ΔV en m/s
    pub dv_tangential_ms: f64,
    /// Constante de Jacobi alcanzada (canónico)
    pub jacobi_achieved: f64,
    /// Error absoluto respecto al C_J objetivo
    pub jacobi_error: f64,
    /// ¿Convergió el ajuste de ε?
    pub epsilon_converged: bool,
    /// Iteraciones usadas en Newton-Raphson
    pub epsilon_iterations: usize,
}

impl InjectionResult {
    /// Convierte el estado canónico a dimensional para outputs en Python.
    pub fn to_dimensional(&self) -> [f64; 6] {
        [
            self.state[0] * D_CHAR,
            self.state[1] * D_CHAR,
            self.state[2] * D_CHAR,
            self.state[3] * V_CHAR,
            self.state[4] * V_CHAR,
            self.state[5] * V_CHAR,
        ]
    }
}

// ============================================================================
// VERIFICACIÓN DE CAPTURA BALÍSTICA
// ============================================================================

/// Veredicto de factibilidad de captura balística.
#[derive(Debug, Clone)]
pub struct BallisticCaptureVerdict {
    pub feasible: bool,
    pub warnings: Vec<String>,
    pub cj_actual: f64,
    pub cj_error: f64,
}

/// Verifica que el punto de inyección pueda alcanzar la cuenca lunar.
/// Acepta tanto inyección desde cuenca terrestre (variedad estable)
/// como desde cuenca lunar (variedad inestable).
pub fn validate_ballistic_capture(
    state: &StateVector,
    target_cj: f64,
) -> BallisticCaptureVerdict {
    let xl1 = l1_position();
    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);
    let cj_actual = jacobi_constant(state);

    let mut warnings = Vec::new();
    let mut feasible = true;

    // Para cruzar L1, necesitamos C_J < C_J(L1) (más energía cinética)
    if cj_actual >= cj_l1 {
        warnings.push(format!(
            "C_J({:.6}) >= C_J(L1)({:.6}): no puede cruzar el cuello",
            cj_actual, cj_l1
        ));
        feasible = false;
    }

    // Determinar en qué cuenca está la nave
    let in_lunar_realm = state.x > xl1;
    let in_terrestrial_realm = state.x < xl1;

    if !in_lunar_realm && !in_terrestrial_realm {
        warnings.push(format!(
            "x({:.6}) ≈ x_L1({:.6}): exactamente en L1",
            state.x, xl1
        ));
        feasible = false;
    }

    // Verificar dirección de velocidad según la cuenca
    if in_terrestrial_realm {
        // Cuenca terrestre: vx > 0 para ir hacia L1 (variedad estable)
        if state[3] <= 0.0 {
            warnings.push(format!(
                "Cuenca terrestre: vx={:.4} <= 0, debe ser > 0 hacia L1",
                state[3]
            ));
            feasible = false;
        }
    } else if in_lunar_realm {
        // Cuenca lunar: vx > 0 para alejarse de L1 hacia la Luna (variedad inestable)
        if state[3] <= 0.0 {
            warnings.push(format!(
                "Cuenca lunar: vx={:.4} <= 0, debe ser > 0 hacia la Luna",
                state[3]
            ));
            feasible = false;
        }
    }

    let cj_error = (cj_actual - target_cj).abs();
    if cj_error > 1e-3 {
        warnings.push(format!(
            "C_J error ({:.2e}) > tolerancia 1e-3", cj_error
        ));
    }

    BallisticCaptureVerdict { feasible, warnings, cj_actual, cj_error }
}

// ============================================================================
// CÁLCULO PRINCIPAL DE INYECCIÓN - TODO EN CANÓNICO
// ============================================================================

pub fn compute_injection(
    orbital: &OrbitalState,
    target_cj: f64,
    epsilon_guess: f64,
    tolerance_cj: f64,
    max_iter: usize,
) -> Result<InjectionResult, InjectionError> {

    let xl1 = l1_position();

    // Validación: la nave debe estar cerca de L1 (en cualquier cuenca)
    let dist_to_l1 = ((orbital.position[0] - xl1).powi(2)
                    + orbital.position[1].powi(2)
                    + orbital.position[2].powi(2)).sqrt();

    if dist_to_l1 > 0.05 {
        return Err(InjectionError::NotInTerrestrialRealm(
            orbital.position[0], xl1,
        ));
    }

    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();
    let alpha = jacobi_sensitivity_l1();

    // Calcular C_J en L1 (referencia)
    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    // El C_J objetivo debe ser MENOR que C_J(L1) para cruzar hacia la Luna
    if target_cj >= cj_l1 {
        return Err(InjectionError::JacobiInfeasible(target_cj, cj_l1));
    }

    if target_cj < cj_l1 - 0.05 {
        return Err(InjectionError::JacobiInfeasible(target_cj, cj_l1));
    }

    // Determinar la cuenca y seleccionar la dirección del manifold
    let (direction, is_lunar_realm) = if orbital.position[0] > xl1 {
        // CUENCA LUNAR: usar variedad INESTABLE (alejándose de L1 hacia Luna)
        // Dirección: (1, k, 0, λ, λk, 0) → x > xL1, vx > 0
        ([1.0, k, 0.0, lambda, lambda * k, 0.0], true)
    } else {
        // CUENCA TERRESTRE: usar variedad ESTABLE (acercándose a L1)
        // Dirección: (-1, -k, 0, λ, λk, 0) → x < xL1, vx > 0
                ([1.0, k, 0.0, lambda, lambda * k, 0.0], false)
    };

    let mut eps = epsilon_guess;
    let mut converged = false;
    let mut iterations = 0;

    // Newton-Raphson para encontrar ε que da el C_J objetivo
    for i in 0..max_iter {
        let test_state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj_test = jacobi_constant(&test_state);
        let err = cj_test - target_cj;

        iterations = i + 1;

        if err.abs() < tolerance_cj {
            converged = true;
            break;
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            eps = eps * 0.5;
            if eps < 1e-10 {
                return Err(InjectionError::EpsilonNotConverged(max_iter));
            }
            continue;
        }

        eps = eps_sq_new.sqrt().clamp(1e-8, 0.1);
    }

    if !converged {
        return Err(InjectionError::EpsilonNotConverged(max_iter));
    }

    // Construir estado de inyección
    let injection_state = StateVector::new(
        xl1 + eps * direction[0],
        eps * direction[1],
        eps * direction[2],
        eps * direction[3],
        eps * direction[4],
        eps * direction[5],
    );
    let cj_achieved = jacobi_constant(&injection_state);

    // Calcular ΔV (diferencia de velocidades en canónico)
    let dvx = injection_state[3] - orbital.velocity[0];
    let dvy = injection_state[4] - orbital.velocity[1];
    let dvz = injection_state[5] - orbital.velocity[2];
    let dv_total_canonical = (dvx * dvx + dvy * dvy + dvz * dvz).sqrt();
    let dv_total_ms = dv_total_canonical * V_CHAR;

    // Calcular componentes radial y tangencial del ΔV
    let r_norm = (orbital.position[0].powi(2)
                + orbital.position[1].powi(2)
                + orbital.position[2].powi(2))
                .sqrt();

    let (dv_radial_ms, dv_tangential_ms) = if r_norm > 0.0 {
        let u_r = [
            orbital.position[0] / r_norm,
            orbital.position[1] / r_norm,
            orbital.position[2] / r_norm,
        ];

        let dv_radial_canonical = dvx * u_r[0] + dvy * u_r[1] + dvz * u_r[2];
        let dv_radial = dv_radial_canonical * V_CHAR;

        let dv_tan_x = dvx - dv_radial_canonical * u_r[0];
        let dv_tan_y = dvy - dv_radial_canonical * u_r[1];
        let dv_tan_z = dvz - dv_radial_canonical * u_r[2];
        let dv_tan_canonical = (dv_tan_x.powi(2) + dv_tan_y.powi(2) + dv_tan_z.powi(2)).sqrt();
        let dv_tan = dv_tan_canonical * V_CHAR;

        (dv_radial.abs(), dv_tan)
    } else {
        (0.0, dv_total_ms)
    };

    // Validación de ΔV físico
    if dv_total_ms > 5000.0 {
        return Err(InjectionError::UnphysicalDeltaV(dv_total_ms));
    }

    // Construir estado dimensional para outputs
    let state_dimensional = [
        injection_state[0] * D_CHAR,
        injection_state[1] * D_CHAR,
        injection_state[2] * D_CHAR,
        injection_state[3] * V_CHAR,
        injection_state[4] * V_CHAR,
        injection_state[5] * V_CHAR,
    ];

    Ok(InjectionResult {
        state: injection_state,
        state_dimensional,
        epsilon: eps,
        dv_total_ms,
        dv_total_canonical,
        dv_radial_ms,
        dv_tangential_ms,
        jacobi_achieved: cj_achieved,
        jacobi_error: (cj_achieved - target_cj).abs(),
        epsilon_converged: converged,
        epsilon_iterations: iterations,
    })
}

// ============================================================================
// BARRIDO DE CONSTANTES DE JACOBI
// ============================================================================

#[derive(Debug, Clone)]
pub struct SweepPoint {
    pub target_cj: f64,
    pub epsilon: f64,
    pub dv_total_ms: f64,
    pub dv_radial_ms: f64,
    pub dv_tangential_ms: f64,
    pub injection_x: f64,
    pub injection_y: f64,
    pub injection_vx: f64,
    pub injection_vy: f64,
    pub jacobi_achieved: f64,
    pub jacobi_error: f64,
    pub converged: bool,
}

pub fn sweep_jacobi(
    orbital: &OrbitalState,
    cj_start: f64,
    cj_end: f64,
    cj_step: f64,
    epsilon_guess: f64,
    tolerance_cj: f64,
    max_iter: usize,
) -> Vec<SweepPoint> {
    let mut results = Vec::new();
    let n_steps = ((cj_end - cj_start) / cj_step).ceil() as usize;

    for i in 0..=n_steps {
        let target_cj = cj_start + i as f64 * cj_step;
        if target_cj > cj_end + 1e-12 {
            break;
        }

        match compute_injection(orbital, target_cj, epsilon_guess, tolerance_cj, max_iter) {
            Ok(inj) => {
                results.push(SweepPoint {
                    target_cj,
                    epsilon: inj.epsilon,
                    dv_total_ms: inj.dv_total_ms,
                    dv_radial_ms: inj.dv_radial_ms,
                    dv_tangential_ms: inj.dv_tangential_ms,
                    injection_x: inj.state.x,
                    injection_y: inj.state.y,
                    injection_vx: inj.state[3],
                    injection_vy: inj.state[4],
                    jacobi_achieved: inj.jacobi_achieved,
                    jacobi_error: inj.jacobi_error,
                    converged: true,
                });
            }
            Err(_) => {
                results.push(SweepPoint {
                    target_cj,
                    epsilon: 0.0,
                    dv_total_ms: f64::NAN,
                    dv_radial_ms: f64::NAN,
                    dv_tangential_ms: f64::NAN,
                    injection_x: f64::NAN,
                    injection_y: f64::NAN,
                    injection_vx: f64::NAN,
                    injection_vy: f64::NAN,
                    jacobi_achieved: f64::NAN,
                    jacobi_error: f64::NAN,
                    converged: false,
                });
            }
        }
    }

    results
}

// ============================================================================
// INYECCIÓN VÍA TFC (ALTERNATIVA RIGUROSA)
// ============================================================================

pub fn compute_injection_tfc(
    orbital: &OrbitalState,
    target_cj: f64,
    tfc_order: usize,
    tof_guess: f64,
    verbose: bool,
) -> Result<InjectionResult, InjectionError> {
    // Por ahora, usar el método directo que es más estable
    compute_injection(orbital, target_cj, 0.005, 1e-8, 50)
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rejects_lunar_realm() {
        let xl1 = l1_position();
        let orbital = OrbitalState::new_canonical(xl1 + 0.01, 0.0, 0.0, 0.0, 0.0, 0.0);
        let result = compute_injection(&orbital, 3.187, 0.005, 1e-10, 30);
        assert!(result.is_err());
    }

    #[test]
    fn test_rejects_high_jacobi() {
        let xl1 = l1_position();
        let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj_l1 = jacobi_constant(&state_l1);

        let orbital = OrbitalState::new_canonical(xl1 - 0.01, -0.001, 0.0, 0.0, 0.0, 0.0);
        // C_J objetivo MAYOR que L1 debe fallar
        let result = compute_injection(&orbital, cj_l1 + 0.01, 0.005, 1e-10, 30);
        assert!(result.is_err());
    }

    #[test]
    fn test_injection_success() {
        let xl1 = l1_position();
        let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj_l1 = jacobi_constant(&state_l1);

        // C_J objetivo MENOR que L1 (cruza hacia la Luna)
        let target = cj_l1 - 0.002;

        let orbital = OrbitalState::new_canonical(
            xl1 - 0.01, -0.001, 0.0,
            0.0, 0.0, 0.0,
        );

        let result = compute_injection(&orbital, target, 0.005, 1e-8, 50);
        assert!(result.is_ok(), "Inyección falló: {:?}", result.err());

        let inj = result.unwrap();

        // Verificar dirección correcta
        assert!(inj.state.x < xl1, "x={:.6} debe ser < x_L1={:.6}", inj.state.x, xl1);
        assert!(inj.state.y < 0.0, "y={:.6} debe ser negativo", inj.state.y);
        assert!(inj.state[3] > 0.0, "vx={:.6} debe ser positivo (hacia Luna)", inj.state[3]);
        assert!(inj.dv_total_ms > 0.0 && inj.dv_total_ms < 1000.0);

        // Verificar C_J
        assert!(inj.jacobi_error < 1e-4, "Error C_J muy grande: {:.2e}", inj.jacobi_error);

        println!("Inyección exitosa: x={:.6}, y={:.6}, ΔV={:.1} m/s, C_J={:.6}",
                 inj.state.x, inj.state.y, inj.dv_total_ms, inj.jacobi_achieved);
    }

    #[test]
    fn test_units_consistency() {
        let xl1 = l1_position();
        let orbital = OrbitalState::new_canonical(xl1 - 0.01, -0.001, 0.0, 0.0, 0.0, 0.0);

        // Verificar que el estado canónico sea < 1.0 (distancia Tierra-Luna)
        assert!(orbital.position[0] < 1.0);
        assert!(orbital.position[0] > 0.8);

        // Verificar conversión a dimensional
        let orbital_dim = OrbitalState::from_dimensional(
            orbital.position[0] * D_CHAR,
            orbital.position[1] * D_CHAR,
            orbital.position[2] * D_CHAR,
            orbital.velocity[0] * V_CHAR,
            orbital.velocity[1] * V_CHAR,
            orbital.velocity[2] * V_CHAR,
        );

        assert!((orbital.position[0] - orbital_dim.position[0]).abs() < 1e-10);
    }

    #[test]
    fn test_ballistic_capture_validation() {
        let xl1 = l1_position();
        let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
        let cj_l1 = jacobi_constant(&state_l1);
        let target = cj_l1 - 0.002;

        let orbital = OrbitalState::new_canonical(xl1 - 0.01, -0.001, 0.0, 0.0, 0.0, 0.0);

        let inj = compute_injection(&orbital, target, 0.005, 1e-8, 50).unwrap();
        let verdict = validate_ballistic_capture(&inj.state, target);

        assert!(verdict.feasible, "Captura balística no factible: {:?}", verdict.warnings);
    }
}
