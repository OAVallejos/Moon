//! propulsion.rs Modelo de propulsión iónica para ASTRO-TFC
//!
//! Implementa la ecuación del cohete de Tsiolkovski para motores iónicos
//! de bajo empuje, con márgenes de ingeniería configurables.
//!
//! PARÁMETROS NOMINALES (configurables):
//! - Empuje: 0.08 N
//! - I_sp: 1,500 s
//! - v_e: I_sp · g₀ = 14,709.975 m/s
//! - Potencia: 1.9 kW
//! - Propelente: Xenón
//!
//! CONTRATO CON lunar.rs / return.rs:
//! Recibe `CaptureDetection` o `EarthCapture` y calcula:
//! - Masa de propelente requerida
//! - Tiempo de encendido
//! - Masa final después de la maniobra

// ============================================================================
// CONFIGURACIÓN DEL MOTOR
// ============================================================================

/// Configuración del motor iónico.
#[derive(Clone)]
pub struct IonEngineConfig {
    /// Empuje nominal [N]
    pub thrust: f64,
    /// Impulso específico [s]
    pub specific_impulse: f64,
    /// Aceleración gravitatoria estándar [m/s²]
    pub g0: f64,
    /// Masa inicial de la nave [kg]
    pub initial_mass_kg: f64,
    /// Margen de ingeniería sobre ΔV [fracción, ej: 0.20 = 20%]
    pub dv_margin: f64,
    /// Margen de ingeniería sobre propelente [fracción]
    pub propellant_margin: f64,
    /// ¿Modo verboso?
    pub verbose: bool,
}

impl Default for IonEngineConfig {
    fn default() -> Self {
        IonEngineConfig {
            thrust: 0.08,
            specific_impulse: 1500.0,
            g0: 9.80665,
            initial_mass_kg: 250.0,
            dv_margin: 0.20,
            propellant_margin: 0.15,
            verbose: false,
        }
    }
}

impl IonEngineConfig {
    /// Velocidad de escape del motor [m/s].
    pub fn exhaust_velocity(&self) -> f64 {
        self.specific_impulse * self.g0
    }
}

// ============================================================================
// RESULTADOS DE MANIOBRA
// ============================================================================

/// Resultado de una maniobra de propulsión.
#[derive(Debug, Clone)]
pub struct ManeuverResult {
    /// ΔV nominal de la maniobra [m/s]
    pub dv_nominal_ms: f64,
    /// ΔV con margen de ingeniería [m/s]
    pub dv_with_margin_ms: f64,
    /// Masa de propelente (mínimo físico) [kg]
    pub propellant_min_kg: f64,
    /// Masa de propelente con márgenes [kg]
    pub propellant_with_margins_kg: f64,
    /// Masa final después de la maniobra [kg]
    pub final_mass_kg: f64,
    /// Tiempo de encendido del motor [días]
    pub burn_time_days: f64,
    /// Fracción de masa consumida [%]
    pub mass_fraction_pct: f64,
    /// ¿La maniobra es factible con el propelente disponible?
    pub is_feasible: bool,
}

// ============================================================================
// ECUACIÓN DE TSIOLKOVSKI
// ============================================================================

/// Calcula la masa final según la ecuación del cohete.
///
/// m_f = m₀ · exp(-ΔV / v_e)
pub fn tsiolkovsky_final_mass(initial_mass_kg: f64, dv_ms: f64, exhaust_velocity_ms: f64) -> f64 {
    initial_mass_kg * (-dv_ms / exhaust_velocity_ms).exp()
}

/// Calcula la masa de propelente requerida.
///
/// m_p = m₀ · (1 - exp(-ΔV / v_e))
pub fn propellant_mass(initial_mass_kg: f64, dv_ms: f64, exhaust_velocity_ms: f64) -> f64 {
    initial_mass_kg - tsiolkovsky_final_mass(initial_mass_kg, dv_ms, exhaust_velocity_ms)
}

/// Calcula el ΔV a partir de la masa de propelente.
///
/// ΔV = v_e · ln(m₀ / m_f)
pub fn delta_v_from_mass(initial_mass_kg: f64, final_mass_kg: f64, exhaust_velocity_ms: f64) -> f64 {
    exhaust_velocity_ms * (initial_mass_kg / final_mass_kg).ln()
}

// ============================================================================
// CÁLCULO DE MANIOBRA
// ============================================================================

/// Calcula los parámetros de una maniobra de propulsión iónica.
///
/// # Argumentos
/// * `dv_required_ms` — ΔV requerido por la maniobra [m/s].
/// * `current_mass_kg` — Masa actual de la nave [kg].
/// * `config` — Configuración del motor.
///
/// # Retorna
/// * `ManeuverResult` con masas, tiempos y factibilidad.
pub fn compute_maneuver(
    dv_required_ms: f64,
    current_mass_kg: f64,
    config: &IonEngineConfig,
) -> ManeuverResult {
    let ve = config.exhaust_velocity();

    // ΔV con margen
    let dv_with_margin = dv_required_ms * (1.0 + config.dv_margin);

    // Masa de propelente (mínimo físico)
    let propellant_min = propellant_mass(current_mass_kg, dv_required_ms, ve);

    // Masa de propelente con márgenes
    let propellant_margin = propellant_mass(current_mass_kg, dv_with_margin, ve);
    let propellant_total = propellant_margin * (1.0 + config.propellant_margin);

    // Masa final
    let final_mass = current_mass_kg - propellant_total;

    // Tiempo de encendido
    // ΔV = (F/m) · t  →  t = ΔV · m_promedio / F
    let avg_mass = (current_mass_kg + final_mass) / 2.0;
    let burn_time_seconds = dv_with_margin * avg_mass / config.thrust;
    let burn_time_days = burn_time_seconds / 86400.0;

    // Fracción de masa
    let mass_fraction_pct = (propellant_total / current_mass_kg) * 100.0;

    // Factibilidad
    let is_feasible = final_mass > 0.0 && propellant_total < current_mass_kg * 0.95;

    if config.verbose {
        println!("   🚀 Maniobra de propulsión iónica:");
        println!("      ΔV requerido:       {:.2} m/s", dv_required_ms);
        println!("      ΔV con margen:      {:.2} m/s", dv_with_margin);
        println!("      v_e (escape):       {:.1} m/s", ve);
        println!("      Masa inicial:       {:.2} kg", current_mass_kg);
        println!("      Propelente (mín):   {:.3} kg", propellant_min);
        println!("      Propelente (total): {:.3} kg", propellant_total);
        println!("      Masa final:         {:.2} kg", final_mass);
        println!("      Tiempo encendido:   {:.1} días", burn_time_days);
        println!("      Fracción masa:      {:.1}%", mass_fraction_pct);
        println!("      Factible:           {}", if is_feasible { "✅" } else { "❌" });
    }

    ManeuverResult {
        dv_nominal_ms: dv_required_ms,
        dv_with_margin_ms: dv_with_margin,
        propellant_min_kg: propellant_min,
        propellant_with_margins_kg: propellant_total,
        final_mass_kg: final_mass,
        burn_time_days,
        mass_fraction_pct,
        is_feasible,
    }
}

// ============================================================================
// PRESUPUESTO COMPLETO DE MISIÓN
// ============================================================================

/// Presupuesto completo de propelente para la misión.
#[derive(Debug, Clone)]
pub struct MissionBudget {
    /// Masa inicial [kg]
    pub initial_mass_kg: f64,
    /// Etapa 1: Inyección a variedad estable
    pub stage1: ManeuverResult,
    /// Etapa 2: Tránsito balístico (ΔV = 0)
    pub stage2_dv_ms: f64,
    /// Etapa 3: Descenso espiral lunar
    pub stage3: ManeuverResult,
    /// Etapa 4: Escape lunar
    pub stage4: ManeuverResult,
    /// Etapa 5: Tránsito balístico retorno (ΔV = 0)
    pub stage5_dv_ms: f64,
    /// Etapa 6: Circularización terrestre
    pub stage6: ManeuverResult,
    /// Masa final [kg]
    pub final_mass_kg: f64,
    /// Propelente total [kg]
    pub total_propellant_kg: f64,
    /// Tiempo total de misión [días]
    pub total_mission_days: f64,
    /// Eficiencia de masa [%]
    pub mass_efficiency_pct: f64,
}

/// Calcula el presupuesto completo de la misión ASTRO-TFC.
///
/// # Argumentos
/// * `dv_injection_ms` — ΔV de inyección a variedad estable [m/s].
/// * `dv_spiral_ms` — ΔV de descenso espiral lunar [m/s].
/// * `dv_escape_ms` — ΔV de escape lunar [m/s].
/// * `dv_circularization_ms` — ΔV de circularización terrestre [m/s].
/// * `tof_ida_days` — Tiempo de vuelo ida [días].
/// * `tof_return_days` — Tiempo de vuelo retorno [días].
/// * `config` — Configuración del motor.
pub fn compute_mission_budget(
    dv_injection_ms: f64,
    dv_spiral_ms: f64,
    dv_escape_ms: f64,
    dv_circularization_ms: f64,
    tof_ida_days: f64,
    tof_return_days: f64,
    config: &IonEngineConfig,
) -> MissionBudget {
    let mut current_mass = config.initial_mass_kg;

    if config.verbose {
        println!("╔══════════════════════════════════════════╗");
        println!("║   PRESUPUESTO DE MISIÓN ASTRO-TFC       ║");
        println!("╚══════════════════════════════════════════╝");
        println!("   Masa inicial: {:.1} kg", current_mass);
    }

    // Etapa 1: Inyección
    let stage1 = compute_maneuver(dv_injection_ms, current_mass, config);
    current_mass = stage1.final_mass_kg;

    // Etapa 2: Tránsito balístico (ΔV = 0, sin consumo)

    // Etapa 3: Descenso espiral lunar
    let stage3 = compute_maneuver(dv_spiral_ms, current_mass, config);
    current_mass = stage3.final_mass_kg;

    // Etapa 4: Escape lunar
    let stage4 = compute_maneuver(dv_escape_ms, current_mass, config);
    current_mass = stage4.final_mass_kg;

    // Etapa 5: Tránsito balístico retorno (ΔV = 0, sin consumo)

    // Etapa 6: Circularización terrestre
    let stage6 = compute_maneuver(dv_circularization_ms, current_mass, config);
    current_mass = stage6.final_mass_kg;

    let total_propellant = config.initial_mass_kg - current_mass;
    let total_days = tof_ida_days + tof_return_days
                   + stage1.burn_time_days
                   + stage3.burn_time_days
                   + stage4.burn_time_days
                   + stage6.burn_time_days;
    let mass_efficiency = (current_mass / config.initial_mass_kg) * 100.0;

    if config.verbose {
        println!("   ─────────────────────────────────────");
        println!("   Propelente total:    {:.2} kg", total_propellant);
        println!("   Masa final:          {:.2} kg", current_mass);
        println!("   Eficiencia de masa:  {:.1}%", mass_efficiency);
        println!("   Tiempo total misión: {:.1} días", total_days);
    }

    MissionBudget {
        initial_mass_kg: config.initial_mass_kg,
        stage1,
        stage2_dv_ms: 0.0,
        stage3,
        stage4,
        stage5_dv_ms: 0.0,
        stage6,
        final_mass_kg: current_mass,
        total_propellant_kg: total_propellant,
        total_mission_days: total_days,
        mass_efficiency_pct: mass_efficiency,
    }
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tsiolkovsky_equation() {
        let m0 = 250.0;
        let dv = 200.0;
        let ve = 1500.0 * 9.80665;

        let mf = tsiolkovsky_final_mass(m0, dv, ve);
        let mp = propellant_mass(m0, dv, ve);

        assert!(mf < m0);
        assert!(mp > 0.0);
        assert!((mf + mp - m0).abs() < 1e-6);

        // Verificar inversa
        let dv_back = delta_v_from_mass(m0, mf, ve);
        assert!((dv_back - dv).abs() < 1e-6);
    }

    #[test]
    fn test_maneuver_feasibility() {
        let config = IonEngineConfig::default();
        let result = compute_maneuver(200.0, 250.0, &config);

        println!("Propulsión: ΔV=200 m/s");
        println!("  Propelente: {:.3} kg", result.propellant_with_margins_kg);
        println!("  Masa final: {:.2} kg", result.final_mass_kg);

        assert!(result.is_feasible);
        assert!(result.propellant_with_margins_kg < 10.0);
        assert!(result.final_mass_kg > 240.0);
    }

    #[test]
    fn test_mission_budget_nominal() {
        let config = IonEngineConfig {
            verbose: true,
            ..IonEngineConfig::default()
        };

        let budget = compute_mission_budget(
            200.87,  // Inyección
            197.7,   // Espiral lunar
            647.83,  // Escape lunar
            0.2,     // Circularización terrestre
            6.7,     // TOF ida
            6.8,     // TOF retorno
            &config,
        );

        println!("\nPresupuesto nominal ASTRO-TFC:");
        println!("  Propelente total: {:.2} kg", budget.total_propellant_kg);
        println!("  Masa final:       {:.2} kg", budget.final_mass_kg);
        println!("  Eficiencia:       {:.1}%", budget.mass_efficiency_pct);
        println!("  Tiempo total:     {:.1} días", budget.total_mission_days);

        assert!(budget.total_propellant_kg < 30.0);
        assert!(budget.final_mass_kg > 220.0);
        assert!(budget.mass_efficiency_pct > 85.0);
    }
}
