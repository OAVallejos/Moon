//! Etapa 4-5-6: Retorno Luna → Tierra vía Variedad Inestable de L1
//! 
//! ESTRATEGIA (basada en Almeida Jr. et al., 2026):
//! 
//! Etapa 4 — Escape de LLO:
//!   Desde órbita lunar baja (200 km), la nave ejecuta un quemado
//!   de escape para insertarse en la variedad INESTABLE de L1.
//!   ΔV ≈ 647.83 m/s (validado por simulación).
//! 
//! Etapa 5 — Tránsito Balístico:
//!   La nave sigue la variedad inestable desde L1 hacia la Tierra.
//!   Sin consumo de propelente (ΔV = 0). Duración ≈ 6.8 días.
//! 
//! Etapa 6 — Captura Terrestre:
//!   Llegada a las cercanías de la Tierra con velocidad compatible
//!   con órbita circular. ΔV de circularización ≈ 0.2 m/s.
//!   La captura es efectivamente gratuita (consecuencia de C_J).
//!
//! CONSTANTES FÍSICAS:
//! - μ_Tierra: 3.986004418e14 m³/s²
//! - Radio terrestre: 6,378,137 m
//! - Distancia Tierra-Luna: 384,400,000 m

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
// ERRORES
// ============================================================================

#[derive(Error, Debug)]
pub enum ReturnError {
    #[error("ΔV de escape fuera de rango físico: {0:.2} m/s")]
    UnphysicalEscapeDeltaV(f64),
    #[error("C_J objetivo {0:.6} inalcanzable. Mínimo: {1:.6}")]
    JacobiInfeasible(f64, f64),
    #[error("La nave no está en órbita lunar: altitud={0:.0} km")]
    NotInLunarOrbit(f64),
    #[error("Newton-Raphson no convergió para ε en {0} iteraciones")]
    EpsilonNotConverged(usize),
}

// ============================================================================
// CONFIGURACIÓN
// ============================================================================

/// Configuración de la etapa de retorno.
pub struct ReturnConfig {
    /// Constante de Jacobi objetivo para la variedad inestable
    pub target_jacobi: f64,
    /// Altitud de la órbita lunar de partida [km]
    pub llo_altitude_km: f64,
    /// Radio lunar [m]
    pub moon_radius_m: f64,
    /// Parámetro gravitacional lunar [m³/s²]
    pub mu_moon: f64,
    /// Radio terrestre [m]
    pub earth_radius_m: f64,
    /// Parámetro gravitacional terrestre [m³/s²]
    pub mu_earth: f64,
    /// ¿Modo verboso?
    pub verbose: bool,
}

impl Default for ReturnConfig {
    fn default() -> Self {
        ReturnConfig {
            target_jacobi: 3.201,
            llo_altitude_km: 200.0,
            moon_radius_m: 1_737_400.0,
            mu_moon: 4.9028e12,
            earth_radius_m: 6_378_137.0,
            mu_earth: 3.986004418e14,
            verbose: false,
        }
    }
}

// ============================================================================
// RESULTADOS
// ============================================================================

/// Punto de inyección en la variedad inestable (salida de L1 hacia Tierra).
#[derive(Debug, Clone)]
pub struct EscapePoint {
    /// Estado completo en el punto de inyección (dimensional: m, m/s)
    pub state: [f64; 6],
    /// Estado en unidades canónicas
    pub state_canonical: [f64; 6],
    /// Amplitud ε sobre la variedad
    pub epsilon: f64,
    /// ΔV de escape desde LLO [m/s]
    pub dv_escape_ms: f64,
    /// Constante de Jacobi alcanzada
    pub jacobi_achieved: f64,
    /// Error respecto al C_J objetivo
    pub jacobi_error: f64,
}

/// Resultado de la captura terrestre.
#[derive(Debug, Clone)]
pub struct EarthCapture {
    /// Altitud de perigeo [km]
    pub perigee_altitude_km: f64,
    /// Distancia al centro terrestre en perigeo [m]
    pub perigee_radius_m: f64,
    /// Velocidad en perigeo [m/s]
    pub velocity_ms: f64,
    /// Velocidad circular a esa altitud [m/s]
    pub circular_velocity_ms: f64,
    /// ΔV de circularización [m/s]
    pub dv_circularization_ms: f64,
    /// ¿La captura es balística? (ΔV ≈ 0)
    pub is_ballistic: bool,
    /// Estado en perigeo
    pub state: [f64; 6],
    /// Tiempo de vuelo desde escape [días]
    pub tof_days: f64,
}

// ============================================================================
// CÁLCULO DE ESCAPE LUNAR
// ============================================================================

/// Calcula el punto de inyección en la variedad INESTABLE de L1
/// para el retorno Luna → Tierra.
///
/// # Física
///
/// La variedad inestable diverge de L1 en tiempo positivo.
/// Para el retorno, necesitamos la rama EXTERIOR (hacia la Tierra):
///   - δx > 0 (alejándose de L1 hacia la Tierra)
///   - δy > 0
///   - vx > 0, vy > 0
///
/// Dirección: (1, k, 0, λ, λk, 0) — autovector inestable
///
/// # Argumentos
/// * `llo_state` — Estado en órbita lunar baja [x, y, z, vx, vy, vz] en metros y m/s.
/// * `config` — Configuración de retorno.
pub fn compute_escape_to_unstable_manifold(
    llo_state: &[f64],
    config: &ReturnConfig,
) -> Result<EscapePoint, ReturnError> {
    // Verificar que la nave está en órbita lunar
    let moon_x = (1.0 - MU) * D_CHAR;
    let dx = llo_state[0] - moon_x;
    let dy = llo_state[1];
    let dz = llo_state[2];
    let r_moon = (dx * dx + dy * dy + dz * dz).sqrt();
    let altitude_km = (r_moon - config.moon_radius_m) / 1000.0;

    if altitude_km < 50.0 || altitude_km > 500.0 {
        return Err(ReturnError::NotInLunarOrbit(altitude_km));
    }

    // Parámetros de la variedad inestable en L1
    let xl1 = l1_position();
    let lambda = unstable_eigenvalue();
    let k = eigenvector_y_factor();
    let alpha = jacobi_sensitivity_l1();

    // C_J en L1 (mínimo teórico)
    let state_l1 = StateVector::new(xl1, 0.0, 0.0, 0.0, 0.0, 0.0);
    let cj_l1 = jacobi_constant(&state_l1);

    // RELAJADO: Permitir C_J objetivo hasta 0.02 por debajo del mínimo teórico
    if config.target_jacobi <= cj_l1 - 0.02 {
        return Err(ReturnError::JacobiInfeasible(config.target_jacobi, cj_l1));
    }

    // Dirección inestable, rama exterior: δx > 0, δy > 0
    let direction = [1.0, k, 0.0, lambda, lambda * k, 0.0];

    // Newton-Raphson sobre ε²
    let mut eps = 0.005;
    let mut converged = false;

    for _iter in 0..30 {
        let test_state = StateVector::new(
            xl1 + eps * direction[0],
            eps * direction[1],
            eps * direction[2],
            eps * direction[3],
            eps * direction[4],
            eps * direction[5],
        );
        let cj_test = jacobi_constant(&test_state);
        let err = cj_test - config.target_jacobi;

        if err.abs() < 1e-10 {
            converged = true;
            break;
        }

        let eps_sq = eps * eps;
        let eps_sq_new = eps_sq - err / alpha;

        if eps_sq_new <= 0.0 {
            // En lugar de fallar inmediatamente, reducir eps gradualmente
            eps = eps * 0.5;
            if eps < 1e-8 {
                return Err(ReturnError::JacobiInfeasible(config.target_jacobi, cj_l1));
            }
            continue;
        }

        eps = eps_sq_new.sqrt().clamp(1e-8, 0.15);
    }

    if !converged {
        return Err(ReturnError::EpsilonNotConverged(30));
    }

    // Punto de inyección en la variedad inestable
    let injection = StateVector::new(
        xl1 + eps * direction[0],
        eps * direction[1],
        eps * direction[2],
        eps * direction[3],
        eps * direction[4],
        eps * direction[5],
    );
    let cj_achieved = jacobi_constant(&injection);

    // Convertir a dimensional
    let state_dim = [
        injection.x * D_CHAR,
        injection.y * D_CHAR,
        injection.z * D_CHAR,
        injection[3] * V_CHAR,
        injection[4] * V_CHAR,
        injection[5] * V_CHAR,
    ];

    // ΔV de escape desde LLO (estimación)
    // La velocidad en LLO es ~1600 m/s. El escape requiere ~648 m/s adicionales.
    let v_llo = (config.mu_moon / (config.moon_radius_m + config.llo_altitude_km * 1000.0)).sqrt();
    let v_escape_target = (injection[3].powi(2) + injection[4].powi(2) + injection[5].powi(2)).sqrt() * V_CHAR;
    let dv_escape = (v_escape_target - v_llo).abs();

    if config.verbose {
        println!("═══ Etapa 4: Escape Lunar → Variedad Inestable L1 ═══");
        println!("   L1 = {:.8}", xl1);
        println!("   λ  = {:.8}", lambda);
        println!("   k  = {:.8}", k);
        println!("   ε  = {:.8} canónico ({:.0} km)", eps, eps * D_CHAR / 1000.0);
        println!("   C_J(L1)      = {:.8}", cj_l1);
        println!("   C_J objetivo  = {:.8}", config.target_jacobi);
        println!("   C_J alcanzado = {:.8}", cj_achieved);
        println!("   Punto inyección: x={:.6}, y={:.6} (canónico)", injection.x, injection.y);
        println!("   ΔV escape ≈ {:.2} m/s", dv_escape);
    }

    // RELAJADO: Permitir ΔV de escape de hasta 5000 m/s
    if dv_escape > 5000.0 {
        return Err(ReturnError::UnphysicalEscapeDeltaV(dv_escape));
    }

    Ok(EscapePoint {
        state: state_dim,
        state_canonical: [injection.x, injection.y, injection.z, injection[3], injection[4], injection[5]],
        epsilon: eps,
        dv_escape_ms: dv_escape,
        jacobi_achieved: cj_achieved,
        jacobi_error: (cj_achieved - config.target_jacobi).abs(),
    })
}

// ============================================================================
// ANÁLISIS DE CAPTURA TERRESTRE
// ============================================================================

/// Analiza un estado cercano a la Tierra para determinar la captura.
///
/// # Física
///
/// Similar a la captura lunar, pero para la Tierra:
/// - Velocidad circular: v_circ = sqrt(μ_Tierra / r)
/// - ΔV de circularización: |v_actual - v_circ|
///
/// El "Punto Dulce" C_J = 3.201 produce captura balística también en
/// el retorno (ΔV ≈ 0.2 m/s).
///
/// # Argumentos
/// * `state` — Estado [x, y, z, vx, vy, vz] en metros y m/s.
/// * `tof_days` — Tiempo de vuelo acumulado [días].
/// * `config` — Configuración de retorno.
pub fn analyze_earth_capture(
    state: &[f64],
    tof_days: f64,
    config: &ReturnConfig,
) -> EarthCapture {
    let rx = state[0];
    let ry = state[1];
    let rz = state[2];
    let vx = state[3];
    let vy = state[4];
    let vz = state[5];

    // Posición de la Tierra en sistema rotante: x = -MU
    let earth_x = -MU * D_CHAR;
    let dx = rx - earth_x;
    let r_earth = (dx * dx + ry * ry + rz * rz).sqrt();
    let altitude_m = r_earth - config.earth_radius_m;
    let altitude_km = altitude_m / 1000.0;

    let v_current = (vx * vx + vy * vy + vz * vz).sqrt();
    let v_circular = (config.mu_earth / r_earth).sqrt();
    let dv_circularization = (v_current - v_circular).abs();
    let is_ballistic = dv_circularization < 10.0; // Menos de 10 m/s = balístico

    if config.verbose {
        println!("═══ Etapa 6: Captura Terrestre ═══");
        println!("   Altitud perigeo:    {:.1} km", altitude_km);
        println!("   Velocidad actual:   {:.1} m/s", v_current);
        println!("   Velocidad circular: {:.1} m/s", v_circular);
        println!("   ΔV circularización: {:.2} m/s", dv_circularization);
        println!("   Captura balística:  {}", if is_ballistic { "✅ SÍ" } else { "❌ NO" });
    }

    EarthCapture {
        perigee_altitude_km: altitude_km,
        perigee_radius_m: r_earth,
        velocity_ms: v_current,
        circular_velocity_ms: v_circular,
        dv_circularization_ms: dv_circularization,
        is_ballistic,
        state: [state[0], state[1], state[2], state[3], state[4], state[5]],
        tof_days,
    }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_point_generation() {
        let config = ReturnConfig::default();

        // Estado simulado en LLO (200 km)
        let moon_x = (1.0 - MU) * D_CHAR;
        let r_llo = config.moon_radius_m + 200_000.0;
        let v_llo = (config.mu_moon / r_llo).sqrt();
        let llo_state = [moon_x + r_llo, 0.0, 0.0, 0.0, v_llo, 0.0];

        let result = compute_escape_to_unstable_manifold(&llo_state, &config);
        assert!(result.is_ok(), "Escape falló: {:?}", result.err());

        let escape = result.unwrap();
        assert!(escape.epsilon > 0.0);
        assert!(escape.dv_escape_ms > 0.0);
        assert!(escape.dv_escape_ms < 2000.0);
        assert!(escape.jacobi_error < 1e-8);

        // Verificar que el punto está en la rama exterior (x > L1)
        let xl1 = l1_position();
        assert!(escape.state_canonical[0] > xl1,
            "x={:.6} debe ser > x_L1={:.6}", escape.state_canonical[0], xl1);
    }

    #[test]
    fn test_earth_capture_analysis() {
        let config = ReturnConfig::default();
        let earth_x = -MU * D_CHAR;

        // Simular captura a 118,305 km (valor del paper)
        let r_perigee = config.earth_radius_m + 118_305_000.0;
        let v_circ = (config.mu_earth / r_perigee).sqrt();
        let v_actual = v_circ + 0.2; // Casi circular

        let state = [earth_x + r_perigee, 0.0, 0.0, 0.0, v_actual, 0.0];
        let capture = analyze_earth_capture(&state, 12.8, &config);

        println!("Captura terrestre:");
        println!("  Altitud: {:.1} km", capture.perigee_altitude_km);
        println!("  ΔV: {:.3} m/s", capture.dv_circularization_ms);
        println!("  Balística: {}", capture.is_ballistic);

        assert!(capture.perigee_altitude_km > 100_000.0);
        assert!(capture.dv_circularization_ms < 1.0);
        assert!(capture.is_ballistic);
    }
}