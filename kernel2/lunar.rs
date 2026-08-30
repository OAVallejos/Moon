//! lunar.rs Etapa 2-3: Detección de captura balística lunar
//!
//! Módulo de análisis orbital PURO. Sin lógica de propulsión.
//! Sin modelos de motor iónico. Sin espirales.
//!
//! Responsabilidades:
//! 1. Detectar el perilunio (máximo acercamiento) en una trayectoria CRTBP
//! 2. Calcular la velocidad circular a esa altitud
//! 3. Estimar el ΔV necesario para circularización
//! 4. Validar que la captura es físicamente viable (no impacto, dentro de SOI)
//!
//! El contrato con el módulo de propulsión es `CaptureDetection`.
//!
//! CONSTANTES FÍSICAS:
//! - Radio lunar:     1,737.4 km  (0.004521 canónico)
//! - μ_Luna:          4.9028e12 m³/s²
//! - SOI lunar:       ~66,000 km  (0.17 canónico)

use thiserror::Error;
use crate::constants::{MU_M, R_MOON, D_CHAR, V_CHAR};

// ============================================================================
// ERRORES
// ============================================================================

#[derive(Error, Debug)]
pub enum LunarError {
    #[error("Trayectoria fuera de SOI lunar: altitud={0:.0} km (máx: 200,000 km)")]
    OutsideSOI(f64),
    #[error("Impacto lunar detectado: altitud={0:.1} km (mín: 10 km)")]
    ImpactDetected(f64),
    #[error("Trayectoria no alcanza la Luna: distancia mínima={0:.0} km")]
    NoLunarEncounter(f64),
    #[error("Velocidad de captura excesiva: ΔV={0:.0} m/s (máx: {1:.0} m/s)")]
    ExcessiveDeltaV(f64, f64),
    #[error("La nave no está en órbita lunar: altitud={0:.0} km (esperado: {1:.0} ± {2:.0} km)")]
    NotInLunarOrbit(f64, f64, f64),
}

// ============================================================================
// RESULTADO DE DETECCIÓN DE CAPTURA
// ============================================================================

/// Resultado del análisis de captura balística.
///
/// Este struct es el "contrato" entre la dinámica orbital (CRTBP)
/// y el sistema de propulsión (motor iónico/químico).
#[derive(Debug, Clone)]
pub struct CaptureDetection {
    /// Altitud de perilunio sobre la superficie lunar [km]
    pub altitude_km: f64,
    /// Distancia al centro de la Luna en perilunio [m]
    pub perilune_radius_m: f64,
    /// Velocidad real de la nave en perilunio [m/s]
    pub velocity_ms: f64,
    /// Velocidad circular a esa altitud [m/s]
    pub circular_velocity_ms: f64,
    /// ΔV necesario para circularizar [m/s]
    /// dv = |v_circular - v_actual|
    pub dv_circularization_ms: f64,
    /// ¿La órbita es elíptica (v_actual > v_circular) o hiperbólica?
    pub is_elliptic: bool,
    /// Excentricidad aproximada (0 = circular, <1 = elíptica, >=1 = escape)
    pub eccentricity_estimate: f64,
    /// Constante de Jacobi en el punto de captura (validación de deriva)
    pub jacobi_at_capture: f64,
    /// Estado completo en perilunio [x, y, z, vx, vy, vz] en metros y m/s
    pub state: [f64; 6],
    /// Tiempo de vuelo acumulado hasta la captura [días]
    pub tof_days: f64,
    /// ¿La captura es físicamente viable? (no impacto, ΔV razonable)
    pub is_viable: bool,
    /// Mensaje descriptivo del estado de captura
    pub diagnosis: String,
}

// ============================================================================
// CONFIGURACIÓN DE DETECCIÓN
// ============================================================================

/// Configuración para el análisis de captura.
#[derive(Clone)]
pub struct CaptureConfig {
    /// Radio de la Luna [m]
    pub moon_radius_m: f64,
    /// Parámetro gravitacional lunar [m³/s²]
    pub mu_moon: f64,
    /// Altitud mínima de seguridad [km] (por debajo = impacto)
    pub min_safe_altitude_km: f64,
    /// Altitud máxima para considerar "captura" [km] (por encima = fuera de SOI)
    pub max_capture_altitude_km: f64,
    /// ΔV máximo aceptable para circularización [m/s]
    pub max_dv_circularization_ms: f64,
    /// Tolerancia de altitud para validación de órbita lunar [km] (default: 15.0)
    pub altitude_tolerance_km: f64,
    /// Altitud objetivo para órbita lunar [km] (default: 800.0)
    pub target_lunar_altitude_km: f64,
    /// ¿Modo verboso?
    pub verbose: bool,
}

impl Default for CaptureConfig {
    fn default() -> Self {
        CaptureConfig {
            moon_radius_m: R_MOON,
            mu_moon: MU_M,
            min_safe_altitude_km: 10.0,
            max_capture_altitude_km: 200_000.0,
            max_dv_circularization_ms: 5000.0,
            altitude_tolerance_km: 15.0,      // ← NUEVO: tolerancia ±15 km
            target_lunar_altitude_km: 800.0,  // ← NUEVO: altitud objetivo
            verbose: false,
        }
    }
}

// ============================================================================
// DETECCIÓN DE PERILUNIO EN TRAYECTORIA
// ============================================================================

pub fn find_perilune(trajectory: &[[f64; 6]]) -> Option<(usize, f64, [f64; 6])> {
    if trajectory.is_empty() {
        return None;
    }

    let moon_x_m = (1.0 - 0.01215058560962404) * D_CHAR;

    let mut min_dist = f64::MAX;
    let mut min_idx = 0;

    for (i, state) in trajectory.iter().enumerate() {
        let dx = state[0] - moon_x_m;
        let dy = state[1];
        let dz = state[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        if dist < min_dist {
            min_dist = dist;
            min_idx = i;
        }
    }

    Some((min_idx, min_dist, trajectory[min_idx]))
}

// ============================================================================
// ANÁLISIS DE CAPTURA
// ============================================================================

pub fn analyze_capture(
    state: &[f64],
    jacobi: f64,
    tof_days: f64,
    config: &CaptureConfig,
) -> Result<CaptureDetection, LunarError> {
    let rx = state[0];
    let ry = state[1];
    let rz = state[2];
    let vx = state[3];
    let vy = state[4];
    let vz = state[5];

    let moon_x_m = (1.0 - 0.01215058560962404) * D_CHAR;

    let dx = rx - moon_x_m;
    let r_moon = (dx * dx + ry * ry + rz * rz).sqrt();
    let altitude_m = r_moon - config.moon_radius_m;
    let altitude_km = altitude_m / 1000.0;

    if altitude_km < config.min_safe_altitude_km {
        return Err(LunarError::ImpactDetected(altitude_km));
    }

    if altitude_km > config.max_capture_altitude_km {
        return Err(LunarError::OutsideSOI(altitude_km));
    }

    let v_current = (vx * vx + vy * vy + vz * vz).sqrt();
    let v_circular = (config.mu_moon / r_moon).sqrt();
    let dv_circularization = (v_current - v_circular).abs();
    let is_elliptic = v_current < v_circular * 2.0_f64.sqrt();

    let specific_energy = 0.5 * v_current * v_current - config.mu_moon / r_moon;
    let eccentricity = if specific_energy < 0.0 {
        let angular_momentum = r_moon * v_current;
        let e_sq = 1.0 + 2.0 * specific_energy * angular_momentum * angular_momentum
                   / (config.mu_moon * config.mu_moon);
        if e_sq > 0.0 { e_sq.sqrt() } else { 0.0 }
    } else {
        1.0 + specific_energy.abs() / (config.mu_moon / r_moon)
    };

    if dv_circularization > config.max_dv_circularization_ms {
        return Err(LunarError::ExcessiveDeltaV(
            dv_circularization, config.max_dv_circularization_ms,
        ));
    }

    let is_viable = altitude_km >= config.min_safe_altitude_km
                 && altitude_km <= config.max_capture_altitude_km
                 && dv_circularization <= config.max_dv_circularization_ms
                 && is_elliptic;

    let diagnosis = if is_viable {
        format!("Captura viable: alt={:.0} km, ΔV={:.1} m/s, e={:.4}",
                altitude_km, dv_circularization, eccentricity)
    } else if !is_elliptic {
        format!("Órbita hiperbólica: alt={:.0} km, v={:.1} > v_esc={:.1} m/s",
                altitude_km, v_current, v_circular * 2.0_f64.sqrt())
    } else {
        format!("Captura no viable: alt={:.0} km, ΔV={:.1} m/s",
                altitude_km, dv_circularization)
    };

    if config.verbose {
        println!("   🌙 Análisis de captura lunar:");
        println!("      Altitud:           {:.1} km", altitude_km);
        println!("      Distancia centro:  {:.1} km", r_moon / 1000.0);
        println!("      Velocidad actual:  {:.1} m/s", v_current);
        println!("      Velocidad circular: {:.1} m/s", v_circular);
        println!("      ΔV circularización: {:.1} m/s", dv_circularization);
        println!("      Excentricidad:     {:.4}", eccentricity);
        println!("      Órbita elíptica:   {}", is_elliptic);
        println!("      Viable:            {}", is_viable);
        println!("      C_J en captura:    {:.8}", jacobi);
        println!("      Diagnóstico:       {}", diagnosis);
    }

    Ok(CaptureDetection {
        altitude_km,
        perilune_radius_m: r_moon,
        velocity_ms: v_current,
        circular_velocity_ms: v_circular,
        dv_circularization_ms: dv_circularization,
        is_elliptic,
        eccentricity_estimate: eccentricity,
        jacobi_at_capture: jacobi,
        state: [rx, ry, rz, vx, vy, vz],
        tof_days,
        is_viable,
        diagnosis,
    })
}

// ============================================================================
// VALIDACIÓN DE ÓRBITA LUNAR (para Stage 4)
// ============================================================================

/// Valida que la nave esté en órbita lunar dentro de la tolerancia
/// 
/// # Argumentos
/// * `state` - Estado dimensional [x, y, z, vx, vy, vz] en metros y m/s
/// * `config` - Configuración de captura
/// 
/// # Retorna
/// * `Ok(())` si la nave está en órbita lunar
/// * `Err(LunarError)` si no está en órbita o está fuera de tolerancia
// Reemplazar estas dos funciones en lunar.rs

pub fn validate_lunar_orbit(
    state: &[f64],
    config: &CaptureConfig,
) -> Result<(), LunarError> {
    let moon_x_m = (1.0 - 0.01215058560962404) * D_CHAR;
    let dx = state[0] - moon_x_m;
    let dy = state[1];
    let dz = state[2];
    let r_moon = (dx * dx + dy * dy + dz * dz).sqrt();
    let altitude_km = (r_moon - config.moon_radius_m) / 1000.0;

    let target = config.target_lunar_altitude_km;
    // Forzar una ventana de tolerancia mínima de 2.0 km si viene más estricta
    let tolerance = config.altitude_tolerance_km.max(2.0); 
    let diff = (altitude_km - target).abs();

    if diff <= tolerance {
        Ok(())
    } else {
        Err(LunarError::NotInLunarOrbit(altitude_km, target, tolerance))
    }
}

pub fn validate_lunar_orbit_with_target(
    state: &[f64],
    target_altitude_km: f64,
    config: &CaptureConfig,
) -> Result<(), LunarError> {
    let moon_x_m = (1.0 - 0.01215058560962404) * D_CHAR;
    let dx = state[0] - moon_x_m;
    let dy = state[1];
    let dz = state[2];
    let r_moon = (dx * dx + dy * dy + dz * dz).sqrt();
    let altitude_km = (r_moon - config.moon_radius_m) / 1000.0;

    let tolerance = config.altitude_tolerance_km.max(2.0);
    let diff = (altitude_km - target_altitude_km).abs();

    if diff <= tolerance {
        Ok(())
    } else {
        Err(LunarError::NotInLunarOrbit(altitude_km, target_altitude_km, tolerance))
    }
}

/// Verifica si un estado está en órbita lunar (dentro de tolerancia)
pub fn is_in_lunar_orbit(state: &[f64], config: &CaptureConfig) -> bool {
    validate_lunar_orbit(state, config).is_ok()
}

/// Obtiene la altitud lunar de un estado dimensional
pub fn get_lunar_altitude_km(state: &[f64]) -> f64 {
    let moon_x_m = (1.0 - 0.01215058560962404) * D_CHAR;
    let dx = state[0] - moon_x_m;
    let dy = state[1];
    let dz = state[2];
    let r_moon = (dx * dx + dy * dy + dz * dz).sqrt();
    (r_moon - R_MOON) / 1000.0
}

// ============================================================================
// ANÁLISIS DESDE TRAYECTORIA COMPLETA
// ============================================================================

pub fn detect_capture_from_trajectory(
    trajectory: &[[f64; 6]],
    trajectory_times: &[f64],
    jacobi_initial: f64,
    config: &CaptureConfig,
) -> Result<CaptureDetection, LunarError> {
    let (idx, min_dist, perilune_state) = find_perilune(trajectory)
        .ok_or(LunarError::NoLunarEncounter(0.0))?;

    let altitude_km = (min_dist - config.moon_radius_m) / 1000.0;

    if altitude_km > config.max_capture_altitude_km {
        return Err(LunarError::NoLunarEncounter(altitude_km));
    }

    let tof_days = if idx < trajectory_times.len() {
        trajectory_times[idx]
    } else {
        0.0
    };

    analyze_capture(&perilune_state, jacobi_initial, tof_days, config)
}

// ============================================================================
// UTILIDADES
// ============================================================================

pub fn circular_velocity_at_altitude(altitude_km: f64, mu_moon: f64, moon_radius_m: f64) -> f64 {
    let r = moon_radius_m + altitude_km * 1000.0;
    (mu_moon / r).sqrt()
}

pub fn escape_velocity_at_altitude(altitude_km: f64, mu_moon: f64, moon_radius_m: f64) -> f64 {
    circular_velocity_at_altitude(altitude_km, mu_moon, moon_radius_m) * 2.0_f64.sqrt()
}

/// Convierte altitud en km a distancia al centro lunar en metros
pub fn altitude_to_radius_m(altitude_km: f64) -> f64 {
    R_MOON + altitude_km * 1000.0
}

/// Convierte distancia al centro lunar en metros a altitud en km
pub fn radius_to_altitude_km(radius_m: f64) -> f64 {
    (radius_m - R_MOON) / 1000.0
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular_velocity_llo() {
        let v_circ = circular_velocity_at_altitude(100.0, MU_M, R_MOON);
        assert!(v_circ > 1500.0 && v_circ < 1700.0,
            "v_circ LLO: {:.1} m/s (esperado ~1630)", v_circ);
    }

    #[test]
    fn test_escape_velocity() {
        let v_esc = escape_velocity_at_altitude(0.0, MU_M, R_MOON);
        assert!(v_esc > 2300.0 && v_esc < 2500.0,
            "v_esc superficie: {:.1} m/s (esperado ~2380)", v_esc);
    }

    #[test]
    fn test_capture_from_perilune() {
        let config = CaptureConfig::default();
        let r_peri = R_MOON + 790_000.0;
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let v_actual = (MU_M / r_peri).sqrt() * 1.05;
        let state = [moon_x + r_peri, 0.0, 0.0, 0.0, v_actual, 0.0];
        let result = analyze_capture(&state, 3.18, 6.7, &config);
        assert!(result.is_ok(), "Captura falló: {:?}", result.err());
        let cap = result.unwrap();
        assert!(cap.altitude_km > 100.0 && cap.altitude_km < 1000.0);
        assert!(cap.dv_circularization_ms > 0.0);
        assert!(cap.is_viable);
    }

    #[test]
    fn test_validate_lunar_orbit() {
        let config = CaptureConfig::default();
        let target_alt = config.target_lunar_altitude_km;
        let tolerance = config.altitude_tolerance_km;
        
        // Estado con altitud dentro de tolerancia (800 km)
        let r = R_MOON + target_alt * 1000.0;
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let state = [moon_x + r, 0.0, 0.0, 0.0, 1600.0, 0.0];
        
        let result = validate_lunar_orbit(&state, &config);
        assert!(result.is_ok(), "Validación falló: {:?}", result.err());
        
        // Estado con altitud fuera de tolerancia
        let r_bad = R_MOON + (target_alt + tolerance + 10.0) * 1000.0;
        let state_bad = [moon_x + r_bad, 0.0, 0.0, 0.0, 1600.0, 0.0];
        let result_bad = validate_lunar_orbit(&state_bad, &config);
        assert!(result_bad.is_err());
    }

    #[test]
    fn test_get_lunar_altitude() {
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let alt_km = 800.0;
        let r = R_MOON + alt_km * 1000.0;
        let state = [moon_x + r, 0.0, 0.0, 0.0, 1600.0, 0.0];
        let computed_alt = get_lunar_altitude_km(&state);
        assert!((computed_alt - alt_km).abs() < 0.1);
    }

    #[test]
    fn test_impact_rejected() {
        let config = CaptureConfig::default();
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let state = [moon_x + R_MOON - 1000.0, 0.0, 0.0, 0.0, 1600.0, 0.0];
        let result = analyze_capture(&state, 3.18, 6.7, &config);
        assert!(result.is_err());
    }

    #[test]
    fn test_outside_soi_rejected() {
        let config = CaptureConfig::default();
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let r = R_MOON + 300_000_000.0;
        let state = [moon_x + r, 0.0, 0.0, 0.0, 500.0, 0.0];
        let result = analyze_capture(&state, 3.18, 6.7, &config);
        assert!(result.is_err());
    }

    #[test]
    fn test_find_perilune_in_trajectory() {
        let moon_x = (1.0 - 0.01215058560962404) * D_CHAR;
        let trajectory: Vec<[f64; 6]> = (0..100)
            .map(|i| {
                let frac = i as f64 / 99.0;
                let r = R_MOON + 5000_000.0 - frac * 4000_000.0;
                [moon_x + r, 0.0, 0.0, 0.0, 1400.0 + frac * 200.0, 0.0]
            })
            .collect();
        let result = find_perilune(&trajectory);
        assert!(result.is_some());
    }

    #[test]
    fn test_altitude_conversion() {
        let alt_km = 800.0;
        let r = altitude_to_radius_m(alt_km);
        let alt_back = radius_to_altitude_km(r);
        assert!((alt_back - alt_km).abs() < 0.001);
    }
}
