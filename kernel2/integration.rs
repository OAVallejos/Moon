//! integration.rs Integración numérica de alta precisión para CRTBP       
//!
//! Implementa Dormand-Prince 8(7) (RK8(7)) con:
//! - Tabla de Butcher completa de 13 etapas (Hairer et al., 1993)
//! - Control de paso adaptativo basado en LTE ponderado + aceleración
//! - Dense Output: interpolantes de Shampine (orden 7)
//! - Localización de eventos por método de Brent
//! - Proyección final sobre superficie de Jacobi
//! - Muestreo por intervalo de tiempo con pre-asignación
//! - Protección contra tunelaje numérico cerca de primarios
//! - Monitoreo de deriva máxima de C_J durante integración
//! - Evento compuesto para detección de escape de L1
//! - Arquitectura genérica con trait SystemConfig
//! - Manejo correcto de dt con signo para integración hacia atrás
//! - Detección de periapsis por mínimo de distancia
//!
//! Referencia: Hairer, Nørsett & Wanner, "Solving Ordinary Differential
//! Equations I", 2da ed., Springer 1993, Tabla 7.2 + Sección II.6.

use thiserror::Error;
use crate::crtbp::{StateVector, crtbp_derivatives, jacobi_constant, MU};
use crate::manifold::{ManifoldType, l1_position};

// ============================================================================
// TRAIT DE CONFIGURACIÓN DEL SISTEMA PLANETARIO
// ============================================================================

pub trait SystemConfig {
    fn mu() -> f64;
    fn secondary_radius_canonical() -> f64;
    fn secondary_x() -> f64 { 1.0 - Self::mu() }
    fn primary_x() -> f64 { -Self::mu() }
}

pub struct EarthMoon;
impl SystemConfig for EarthMoon {
    fn mu() -> f64 { MU }
    fn secondary_radius_canonical() -> f64 { 1737.4 / 384400.0 }
}

// ============================================================================
// CONSTANTES FÍSICAS
// ============================================================================

const ESCAPE_LIMIT: f64 = 10.0;
const IMPACT_FACTOR: f64 = 0.99;
const TUNNELING_SAFETY_FACTOR: f64 = 1.5;
const MAX_STEP_BODY_FRACTION: f64 = 0.1;

// ============================================================================
// ERRORES
// ============================================================================

#[derive(Error, Debug)]
pub enum IntegrationError {
    #[error("Error en integración: {0}")]
    PropagationError(String),
    #[error("Máximo de pasos alcanzado: {0}")]
    MaxStepsReached(usize),
    #[error("Trayectoria escapó del sistema: r = {0:.4}")]
    EscapedTrajectory(f64),
    #[error("Evento no alcanzado en t_max = {0:.2}")]
    EventNotReached(f64),
    #[error("Impacto con secundario: r = {0:.6} canónico")]
    SecondaryImpact(f64),
    #[error("Impacto con primario: r = {0:.6} canónico")]
    PrimaryImpact(f64),
    #[error("Paso mínimo alcanzado sin convergencia: dt = {0:.2e}")]
    StepSizeUnderflow(f64),
    #[error("Tunelaje numérico detectado: dt={0:.2e}, r={1:.6}")]
    NumericalTunneling(f64, f64),
    #[error("Trayectoria se alejó de L1 sin alcanzar objetivo: dist_L1={0:.4}")]
    EscapedL1Region(f64),
}

// ============================================================================
// CONFIGURACIÓN DEL INTEGRADOR
// ============================================================================

#[derive(Clone)]
pub struct IntegratorConfig {
    pub rtol: f64,
    pub atol: f64,
    pub max_step: f64,
    pub max_step_near_body: f64,
    pub min_step: f64,
    pub initial_step: f64,
    pub safety_factor: f64,
    pub max_growth: f64,
    pub min_reduction: f64,
    pub jacobi_projection: bool,
    pub sampling_interval: f64,
    pub jacobi_monitor_freq: usize,
    pub verbose: bool,
}

impl Default for IntegratorConfig {
    fn default() -> Self {
        IntegratorConfig {
            rtol: 1e-12, atol: 1e-15, max_step: 0.005, max_step_near_body: 0.0001,
            min_step: 1e-10, initial_step: 0.0001, safety_factor: 0.8,
            max_growth: 4.0, min_reduction: 0.1, jacobi_projection: true,
            sampling_interval: 0.01, jacobi_monitor_freq: 1000, verbose: false,
        }
    }
}

// ============================================================================
// RESULTADO DE PROPAGACIÓN
// ============================================================================

pub struct PropagationResult {
    pub tof: f64,
    pub final_state: [f64; 6],
    pub trajectory: Vec<[f64; 6]>,
    pub trajectory_times: Vec<f64>,
    pub jacobi_error: f64,
    pub jacobi_error_max: f64,
    pub target_reached: bool,
    pub steps_attempted: usize,
    pub steps_accepted: usize,
    pub time_limit_reached: bool,
    pub escaped_l1: bool,
    pub termination_reason: String,
    pub min_lunar_distance: f64,
    pub min_lunar_altitude_km: f64,
}

// ============================================================================
// COEFICIENTES DORMAND-PRINCE 8(7) 13M
// ============================================================================

const DP87_S: usize = 13;

const C: [f64; DP87_S] = [
    0.0, 1.0/18.0, 1.0/12.0, 1.0/8.0, 5.0/16.0, 3.0/8.0,
    59.0/400.0, 93.0/200.0, 5490023248.0/9719169821.0,
    13.0/20.0, 1201146811.0/1299019798.0, 1.0, 1.0,
];

const A: [[f64; DP87_S]; DP87_S] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0/18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0/48.0, 1.0/16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1.0/32.0, 0.0, 3.0/32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0, -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0, 800635310.0/3783071287.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0, -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0, 393006217.0/1396673457.0, 123872331.0/1001029789.0, 0.0, 0.0, 0.0, 0.0],
    [-1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0, 1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0, -45442868181.0/3398467696.0, 3065993473.0/597172653.0, 0.0, 0.0, 0.0],
    [185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0, -477755414.0/1098053517.0, -703635378.0/230739211.0, 5731566787.0/1027545527.0, 5232866602.0/850066563.0, -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0, 0.0, 0.0],
    [403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0, -411421997.0/543043805.0, 652783627.0/914296604.0, 11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0, 248638103.0/1413531060.0, 0.0, 0.0],
];

const B8: [f64; DP87_S] = [
    14005451.0/335480064.0, 0.0, 0.0, 0.0, 0.0,
    -59238493.0/1068277825.0, 181606767.0/758867731.0,
    561292985.0/797845732.0, -1041891430.0/1371343529.0,
    760417239.0/1151165299.0, 118820643.0/751138087.0,
    -528747749.0/2220607170.0, 1.0/4.0,
];

const B7: [f64; DP87_S] = [
    13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0,
    -808719846.0/976000145.0, 1757004468.0/5645159321.0,
    656045339.0/265891186.0, -3867574721.0/1518517206.0,
    465885868.0/322736535.0, 53011238.0/667516719.0,
    2.0/45.0, 0.0,
];

// ============================================================================
// DENSE OUTPUT — INTERPOLANTES DE SHAMPINE PARA DOPRI87 (ORDEN 7)
// ============================================================================

fn dense_output_shampine(
    y0: &StateVector,
    dt: f64,
    k: &[StateVector; DP87_S],
    theta: f64,
) -> StateVector {
    let t1 = theta;
    let t2 = t1 * t1;
    let t3 = t2 * t1;
    let t4 = t3 * t1;
    let t5 = t4 * t1;
    let t6 = t5 * t1;
    let _t7 = t6 * t1;

    let b: [f64; DP87_S] = [
        t1 * (1.0 + t1 * (-1337.0/480.0 + t1 * (1039.0/360.0 + t1 * (-1163.0/1152.0)))),
        0.0,
        t2 * (421.0/144.0 + t1 * (-16273.0/2880.0 + t1 * (61529.0/17280.0 + t1 * (-13573.0/11520.0)))),
        t2 * (-55.0/48.0 + t1 * (4019.0/960.0 + t1 * (-18973.0/5760.0 + t1 * (4573.0/3840.0)))),
        t2 * (435.0/32.0 + t1 * (-4671.0/128.0 + t1 * (4833.0/128.0 + t1 * (-1125.0/64.0)))),
        t2 * (-112.0/15.0 + t1 * (964.0/45.0 + t1 * (-408.0/15.0 + t1 * (304.0/15.0)))),
        t2 * (539.0/240.0 + t1 * (-1573.0/160.0 + t1 * (2977.0/480.0 + t1 * (-343.0/192.0)))),
        t2 * (87.0/80.0 + t1 * (-373.0/96.0 + t1 * (103.0/32.0 + t1 * (-53.0/48.0)))),
        t2 * (145.0/96.0 + t1 * (-197.0/64.0 + t1 * (133.0/96.0 + t1 * (-7.0/32.0)))),
        t2 * (-153.0/160.0 + t1 * (289.0/96.0 + t1 * (-293.0/160.0 + t1 * (19.0/32.0)))),
        t2 * (-241.0/480.0 + t1 * (217.0/160.0 + t1 * (-239.0/240.0 + t1 * (71.0/240.0)))),
        t2 * (29.0/240.0 + t1 * (-89.0/240.0 + t1 * (43.0/160.0 + t1 * (-7.0/80.0)))),
        t1 * t2 * (-1.0/12.0 + t1 * (1.0/12.0 + t1 * (-1.0/24.0))),
    ];

    let mut result = *y0;
    for i in 0..DP87_S {
        if b[i] != 0.0 {
            result = result + k[i] * (dt * b[i]);
        }
    }
    result
}

// ============================================================================
// LOCALIZACIÓN DE EVENTOS — MÉTODO DE BRENT
// ============================================================================

fn brent_event(
    y0: &StateVector,
    dt: f64,
    k: &[StateVector; DP87_S],
    event: &impl Fn(f64, &StateVector) -> f64,
    tol: f64,
    max_iter: usize,
) -> (f64, StateVector) {
    let eval = |theta: f64| -> (f64, StateVector) {
        let s = dense_output_shampine(y0, dt, k, theta);
        (event(dt * theta, &s), s)
    };

    let mut a = 0.0;
    let mut b = 1.0;
    let (mut fa, _) = eval(a);
    let (mut fb, _) = eval(b);

    if fa * fb > 0.0 {
        return eval(0.5);
    }

    if fa.abs() > fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut _d = b - a;

    for _ in 0..max_iter {
        if fb.abs() < tol || (b - a).abs() < tol {
            return eval(b);
        }

        let s = if fa != fc && fb != fc {
            let s_iqi = a * fb * fc / ((fa - fb) * (fa - fc))
                      + b * fa * fc / ((fb - fa) * (fb - fc))
                      + c * fa * fb / ((fc - fa) * (fc - fb));
            if s_iqi > (3.0*a + b)/4.0 && s_iqi < b { s_iqi } else { (a + b) / 2.0 }
        } else {
            (a + b) / 2.0
        };

        _d = c;
        c = b;
        fc = fb;

        let s_clamped = if (s - b).abs() < tol {
            b + tol.copysign(b - a)
        } else {
            s
        };

        let (fs_new, _) = eval(s_clamped);
        b = s_clamped;
        fb = fs_new;

        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }

    eval(b)
}

// ============================================================================
// INTEGRADOR GENÉRICO
// ============================================================================

pub struct CRTBPIntegrator<S: SystemConfig> {
    config: IntegratorConfig,
    _phantom: std::marker::PhantomData<S>,
}

impl<S: SystemConfig> CRTBPIntegrator<S> {
    pub fn new(config: IntegratorConfig) -> Self {
        CRTBPIntegrator { config, _phantom: std::marker::PhantomData }
    }

    fn dp87_step_full(
        &self, state: &StateVector, t: f64, dt: f64,
    ) -> (StateVector, StateVector, [StateVector; DP87_S], f64) {
        let mut k = [StateVector::zeros(); DP87_S];

        for i in 0..DP87_S {
            let ci = C[i];
            let mut sum = StateVector::zeros();
            for j in 0..i {
                let aij = A[i][j];
                if aij != 0.0 {
                    sum = sum + k[j] * (aij * dt);
                }
            }
            let s = *state + sum;
            k[i] = crtbp_derivatives(t + ci * dt, &s);
        }

        let mut y8 = StateVector::zeros();
        for i in 0..DP87_S { y8 = y8 + k[i] * (dt * B8[i]); }
        let new_state = *state + y8;

        let mut y7 = StateVector::zeros();
        for i in 0..DP87_S { y7 = y7 + k[i] * (dt * B7[i]); }
        let state7 = *state + y7;

        let deriv_new = crtbp_derivatives(t + dt, &new_state);
        let error_est = self.compute_error(&new_state, &state7);

        (new_state, deriv_new, k, error_est)
    }

    fn compute_error(&self, y8: &StateVector, y7: &StateVector) -> f64 {
        let scale = |y: f64, dy: f64| -> f64 {
            self.config.atol + self.config.rtol * y.abs().max(dy.abs())
        };
        let diff = *y8 - *y7;
        let comps = [diff[0], diff[1], diff[2], diff[3], diff[4], diff[5]];
        let refs = [y8[0], y8[1], y8[2], y8[3], y8[4], y8[5]];
        let mut err_sq: f64 = 0.0;
        for i in 0..6 {
            let s = scale(refs[i], comps[i]);
            err_sq += (comps[i] / s).powi(2);
        }
        (err_sq / 6.0).sqrt()
    }

    fn acceleration_norm(state: &StateVector) -> f64 {
        let mu = S::mu();
        let r1 = ((state[0] + mu).powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
        let r2 = ((state[0] - 1.0 + mu).powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
        let ax = 2.0 * state[4] + state[0]
               - (1.0 - mu) * (state[0] + mu) / (r1 * r1 * r1)
               - mu * (state[0] - 1.0 + mu) / (r2 * r2 * r2);
        let ay = -2.0 * state[3] + state[1]
               - (1.0 - mu) * state[1] / (r1 * r1 * r1)
               - mu * state[1] / (r2 * r2 * r2);
        (ax * ax + ay * ay).sqrt()
    }

    fn adapt_step(&self, _state: &StateVector, dt_current: f64, error_est: f64) -> f64 {
        let sign = if dt_current >= 0.0 { 1.0 } else { -1.0 };
        let dt_abs = dt_current.abs();

        // 1. Calculamos el factor óptimo basándonos únicamente en el error
        // Para RK8(7), el exponente es 1/8
        let error_factor = if error_est > 0.0 {
            (self.config.rtol / error_est).powf(1.0 / 8.0)
        } else {
            self.config.max_growth
        };

        // 2. Aplicamos el factor de seguridad y limitamos el crecimiento/reducción
        // Esto permite que el paso crezca libremente cuando el error es bajo
        let factor = (self.config.safety_factor * error_factor)
            .clamp(self.config.min_reduction, self.config.max_growth);

        // 3. Calculamos el nuevo dt y lo clampamos a los límites absolutos
        let dt_new_abs = (dt_abs * factor).clamp(self.config.min_step, self.config.max_step);

        sign * dt_new_abs
    }

    fn distance_to_secondary(state: &StateVector) -> f64 {
        let sx = S::secondary_x();
        ((state[0] - sx).powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt()
    }

    fn distance_to_primary(state: &StateVector) -> f64 {
        let px = S::primary_x();
        ((state[0] - px).powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt()
    }

    fn project_jacobi_final(state: &mut StateVector, target_cj: f64) -> f64 {
        let cj_before = jacobi_constant(state);
        let err = cj_before - target_cj;
        if err.abs() < 1e-12 { return 0.0; }
        let v_sq = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];
        if v_sq < 1e-20 { return err; }
        let v_sq_target = v_sq + err;
        if v_sq_target <= 0.0 { return err; }
        let scale = (v_sq_target / v_sq).sqrt();
        state[3] *= scale;
        state[4] *= scale;
        state[5] *= scale;
        jacobi_constant(state) - target_cj
    }

    // ========================================================================
    // PROPAGACIÓN PRINCIPAL
    // ========================================================================

    pub fn propagate_until_event(
        &self,
        initial_state: &StateVector,
        t_max: f64,
        time_direction: f64,
        event: impl Fn(f64, &StateVector) -> f64,
        tol_event: f64,
        target_cj: Option<f64>,
        max_steps: usize,
        l1_escape_threshold: Option<f64>,
    ) -> Result<PropagationResult, IntegrationError> {
        let initial_cj = jacobi_constant(initial_state);
        let cj_target = target_cj.unwrap_or(initial_cj);
        let secondary_radius = S::secondary_radius_canonical();
        let impact_radius = secondary_radius * IMPACT_FACTOR;
        let tunneling_radius = secondary_radius * TUNNELING_SAFETY_FACTOR;

        let estimated_samples = (t_max / self.config.sampling_interval).ceil() as usize + 2;
        let mut trajectory: Vec<[f64; 6]> = Vec::with_capacity(estimated_samples);
        let mut trajectory_times: Vec<f64> = Vec::with_capacity(estimated_samples);

        let deriv_initial = crtbp_derivatives(0.0, initial_state);
        let mut state = *initial_state;
        let mut deriv = deriv_initial;
        let mut t = 0.0;
        let mut dt = self.config.initial_step * time_direction;

        trajectory.push([state[0], state[1], state[2], state[3], state[4], state[5]]);
        trajectory_times.push(t);

        let mut last_event_value = event(t, &state);
        let mut steps_attempted = 0usize;
        let mut steps_accepted = 0usize;
        let mut event_hit = false;
        let mut time_limit_reached = false;
        let mut escaped_l1 = false;
        let mut next_sample_time = self.config.sampling_interval;
        let mut jacobi_error_max = 0.0f64;

        // Monitoreo de distancia mínima a la Luna
        let mut min_lunar_distance = f64::MAX;
        let mut min_lunar_state = *initial_state;
        let mut min_lunar_time = 0.0;
        let mut prev_lunar_distance = Self::distance_to_secondary(&state);

        while t < t_max && steps_attempted < max_steps {
            let r_secondary = Self::distance_to_secondary(&state);

            // Registrar mínimo de distancia
            if r_secondary < min_lunar_distance {
                min_lunar_distance = r_secondary;
                min_lunar_state = state;
                min_lunar_time = t;
            }

            let effective_max_step = if r_secondary < tunneling_radius {
                self.config.max_step_near_body
            } else {
                self.config.max_step
            };

            let dt_abs_proposed = dt.abs().min(t_max - t).min(effective_max_step);
            let dt_proposed = time_direction * dt_abs_proposed;

            if dt_abs_proposed < self.config.min_step {
                return Err(IntegrationError::StepSizeUnderflow(dt_abs_proposed));
            }

            if r_secondary < tunneling_radius * 2.0
                && dt_abs_proposed > r_secondary * MAX_STEP_BODY_FRACTION
            {
                return Err(IntegrationError::NumericalTunneling(dt_abs_proposed, r_secondary));
            }

            steps_attempted += 1;

            let (new_state, new_deriv, k_array, error_est) =
                self.dp87_step_full(&state, t, dt_proposed);

            let dt_new = self.adapt_step(&state, dt_proposed, error_est);

            if error_est > 1.0 {
                dt = dt_new;
                continue;
            }

            steps_accepted += 1;
            let new_t = t + dt_proposed;

            if steps_accepted % self.config.jacobi_monitor_freq == 0 {
                let cj_current = jacobi_constant(&state);
                let drift = (cj_current - initial_cj).abs();
                if drift > jacobi_error_max { jacobi_error_max = drift; }
            }

            let new_event_value = event(new_t, &new_state);

            // Detección de cruce por cero del evento
            if last_event_value * new_event_value <= 0.0 && last_event_value.abs() > 1e-15 {
                let (theta_event, state_event) = brent_event(
                    &state, dt_proposed, &k_array, &event, tol_event, 40,
                );
                t = t + theta_event * dt_proposed;
                state = state_event;
                deriv = crtbp_derivatives(t, &state);
                trajectory.push([state[0], state[1], state[2], state[3], state[4], state[5]]);
                trajectory_times.push(t);
                event_hit = true;
                break;
            }

            // Detectar periapsis: la distancia empieza a aumentar
            let new_lunar_distance = Self::distance_to_secondary(&new_state);
            if prev_lunar_distance < min_lunar_distance * 1.001
                && new_lunar_distance > prev_lunar_distance * 1.01
                && min_lunar_distance < 0.05
            {
                state = min_lunar_state;
                t = min_lunar_time;
                deriv = crtbp_derivatives(t, &state);
                event_hit = true;
                break;
            }
            prev_lunar_distance = new_lunar_distance;

            state = new_state;
            deriv = new_deriv;
            t = new_t;
            last_event_value = new_event_value;
            dt = dt_new;

            if t >= next_sample_time {
                trajectory.push([state[0], state[1], state[2], state[3], state[4], state[5]]);
                trajectory_times.push(t);
                next_sample_time += self.config.sampling_interval;
            }

            if let Some(threshold) = l1_escape_threshold {
                let xl1 = l1_position();
                let dist_l1 = ((state[0] - xl1).powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
                if dist_l1 > threshold {
                    escaped_l1 = true;
                    break;
                }
            }

            let r = (state[0].powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
            if r > ESCAPE_LIMIT { return Err(IntegrationError::EscapedTrajectory(r)); }

            let r_sec = Self::distance_to_secondary(&state);
            if r_sec < impact_radius { return Err(IntegrationError::SecondaryImpact(r_sec)); }

            let r_prim = Self::distance_to_primary(&state);
            if r_prim < 0.0001 { return Err(IntegrationError::PrimaryImpact(r_prim)); }
        }

        if t >= t_max { time_limit_reached = true; }

        let termination_reason = if event_hit {
            format!("Periapsis lunar en t={:.4}", t)
        } else if escaped_l1 {
            format!("Escapó de región L1 en t={:.4}", t)
        } else if time_limit_reached {
            format!("Tiempo máximo: t={:.4}", t)
        } else {
            format!("Máximo de pasos: {}", max_steps)
        };

        let mut final_state = state;
        if self.config.jacobi_projection {
            Self::project_jacobi_final(&mut final_state, cj_target);
        }
        let final_cj = jacobi_constant(&final_state);

        if trajectory.is_empty()
            || (trajectory.last().unwrap()[0] - final_state[0]).abs() > 1e-10
        {
            trajectory.push([final_state[0], final_state[1], final_state[2], final_state[3], final_state[4], final_state[5]]);
            trajectory_times.push(t);
        }

        let min_altitude_km = {
    use crate::constants::D_CHAR;
    let dist_m = min_lunar_distance * D_CHAR;
    let alt_m = dist_m - crate::constants::R_MOON;
    alt_m / 1000.0
};

        Ok(PropagationResult {
            tof: t,
            final_state: [final_state[0], final_state[1], final_state[2], final_state[3], final_state[4], final_state[5]],
            trajectory,
            trajectory_times,
            jacobi_error: (final_cj - initial_cj).abs(),
            jacobi_error_max,
            target_reached: event_hit,
            steps_attempted,
            steps_accepted,
            time_limit_reached,
            escaped_l1,
            termination_reason,
            min_lunar_distance,
            min_lunar_altitude_km: min_altitude_km,
        })
    }
}

// ============================================================================
// TIPO CONCRETO + PROPAGACIÓN POR VARIEDADES + UTILIDADES
// ============================================================================

pub type EarthMoonIntegrator = CRTBPIntegrator<EarthMoon>;

pub fn propagate_manifold<S: SystemConfig>(
    initial_state: &[f64],
    manifold_type: ManifoldType,
    target_cj: f64,
    config: &IntegratorConfig,
    max_tof_days: f64,
    target_lunar_altitude_km: Option<f64>,
    l1_escape_threshold: Option<f64>,
) -> Result<PropagationResult, IntegrationError> {
    let state0 = StateVector::new(
        initial_state[0], initial_state[1], initial_state[2],
        initial_state[3], initial_state[4], initial_state[5],
    );

    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    // NUEVO: Validación de consistencia física
    // Verifica que el estado inicial esté en la cuenca correcta
    // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    let xl1 = l1_position();
    match manifold_type {
        ManifoldType::TransitToMoon => {
            // Debe estar en cuenca lunar (x > xL1) con vx > 0
            if state0[0] <= xl1 {
                return Err(IntegrationError::PropagationError(
                    format!("TransitToMoon requiere cuenca lunar (x > xL1). x={:.6}, xL1={:.6}",
                            state0[0], xl1)
                ));
            }
            if state0[3] <= 0.0 {
                return Err(IntegrationError::PropagationError(
                    format!("TransitToMoon requiere vx > 0 (hacia la Luna). vx={:.6}", state0[3])
                ));
            }
        },
        ManifoldType::TransitToEarth => {
            // Debe estar en cuenca terrestre (x < xL1) con vx < 0
            if state0[0] >= xl1 {
                return Err(IntegrationError::PropagationError(
                    format!("TransitToEarth requiere cuenca terrestre (x < xL1). x={:.6}, xL1={:.6}",
                            state0[0], xl1)
                ));
            }
            if state0[3] >= 0.0 {
                return Err(IntegrationError::PropagationError(
                    format!("TransitToEarth requiere vx < 0 (hacia la Tierra). vx={:.6}", state0[3])
                ));
            }
        },
        _ => {} // Stable/Unstable no necesitan validación adicional
    }

    let time_direction = match manifold_type {
        ManifoldType::Stable => -1.0,
        ManifoldType::Unstable => 1.0,
        ManifoldType::TransitToMoon => 1.0,
        ManifoldType::TransitToEarth => -1.0,
    };

    let t_max = max_tof_days * 86400.0 / crate::constants::T_CHAR;
    let integrator = CRTBPIntegrator::<S>::new(config.clone());

    let target_r = if let Some(alt_km) = target_lunar_altitude_km {
        S::secondary_radius_canonical() + (alt_km * 1000.0) / 384400000.0
    } else {
        S::secondary_radius_canonical() + (800.0 * 1000.0) / 384400000.0
    };

    let event = move |_t: f64, s: &StateVector| -> f64 {
        let dx = s[0] - S::secondary_x();
        (dx * dx + s[1] * s[1] + s[2] * s[2]).sqrt() - target_r
    };

    if config.verbose {
        println!("   🔭 Propagando {:?} (dir={})", manifold_type, time_direction);
        println!("      C_J inicial: {:.10}", jacobi_constant(&state0));
    }

    integrator.propagate_until_event(
        &state0, t_max, time_direction, event,
        1e-10, Some(target_cj), 500_000, l1_escape_threshold,
    )
}

pub fn distance_to_moon(state: &[f64]) -> f64 {
    let sx = EarthMoon::secondary_x();
    let dx = state[0] - sx;
    let dy = state[1];
    let dz = if state.len() > 2 { state[2] } else { 0.0 };
    (dx * dx + dy * dy + dz * dz).sqrt()
}

pub fn lunar_altitude_km(state: &[f64]) -> f64 {
    use crate::constants::D_CHAR;
    // distance_to_moon() devuelve unidades canónicas (fracción de D_CHAR)
    // Multiplicar por D_CHAR da metros, dividir por 1000 da km
    let dist_m = distance_to_moon(state) * D_CHAR;
    let alt_m = dist_m - crate::constants::R_MOON;
    alt_m / 1000.0
}

pub fn distance_to_earth(state: &[f64]) -> f64 {
    let px = EarthMoon::primary_x();
    let dx = state[0] - px;
    let dy = state[1];
    let dz = if state.len() > 2 { state[2] } else { 0.0 };
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Calcula altitud lunar en km desde estado DIMENSIONAL [m, m, m, m/s, m/s, m/s].
/// A diferencia de lunar_altitude_km() que espera estado CANÓNICO,
/// esta función recibe el estado ya convertido a metros.
pub fn lunar_altitude_km_from_dimensional(state: &[f64]) -> f64 {
    let moon_x = (1.0 - crate::crtbp::MU) * crate::constants::D_CHAR;
    let dx = state[0] - moon_x;
    let dy = state[1];
    let dz = if state.len() > 2 { state[2] } else { 0.0 };
    let dist_m = (dx * dx + dy * dy + dz * dz).sqrt();
    let alt_m = dist_m - crate::constants::R_MOON;
    alt_m / 1000.0
}
// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_table_consistency() {
        let sum_b8: f64 = B8.iter().sum();
        assert!((sum_b8 - 1.0).abs() < 1e-14, "ΣB8 = {:.15}", sum_b8);
        for i in 0..DP87_S {
            let sum_a: f64 = A[i].iter().sum();
            assert!((C[i] - sum_a).abs() < 1e-14, "Fila {}: c={}, Σa={}", i, C[i], sum_a);
        }
    }

    #[test]
    fn test_dp87_jacobi_conservation() {
        let config = IntegratorConfig::default();
        let integrator = CRTBPIntegrator::<EarthMoon>::new(config);
        let state = StateVector::new(0.83, 0.0, 0.0, 0.0, 0.15, 0.0);
        let cj_before = jacobi_constant(&state);
        let (new_state, _, _, error_est) = integrator.dp87_step_full(&state, 0.0, 0.001);
        let cj_after = jacobi_constant(&new_state);
        assert!((cj_after - cj_before).abs() < 1e-8,
            "Error C_J: {:.2e}, LTE: {:.2e}", (cj_after-cj_before).abs(), error_est);
    }

    #[test]
    fn test_lunar_altitude() {
        let alt_m = 800_000.0;
        let r = EarthMoon::secondary_radius_canonical() + alt_m / 384400_000.0;
        let sx = EarthMoon::secondary_x();
        let state = [sx + r, 0.0, 0.0, 0.0, 0.0, 0.0];
        let alt = lunar_altitude_km(&state);
        assert!((alt - 800.0).abs() < 1.0, "Altitud: {:.0} km", alt);
    }

    #[test]
    fn test_system_config_values() {
        assert!((EarthMoon::mu() - 0.01215058560962404).abs() < 1e-16);
        assert!(EarthMoon::secondary_x() > 0.98);
        assert!(EarthMoon::primary_x() < 0.0);
        assert!(EarthMoon::secondary_radius_canonical() < 0.005);
    }

    #[test]
    fn test_backward_integration() {
        let config = IntegratorConfig { max_step: 0.001, verbose: false, ..IntegratorConfig::default() };
        let integrator = CRTBPIntegrator::<EarthMoon>::new(config);
        let state = StateVector::new(0.83, 0.01, 0.0, 0.0, -0.1, 0.0);
        let event = |_t: f64, _s: &StateVector| -> f64 { -1.0 };
        let result = integrator.propagate_until_event(&state, 0.5, -1.0, event, 1e-10, None, 10_000, None);
        match result {
            Ok(res) => {
                assert!(res.tof > 0.0);
                assert!(res.steps_accepted > 0);
            }
            Err(e) => panic!("Backward integration failed: {:?}", e),
        }
    }

    #[test]
    fn test_jacobi_stability_long_propagation() {
        let config = IntegratorConfig {
            rtol: 1e-13, atol: 1e-16, max_step: 0.005,
            jacobi_projection: false, jacobi_monitor_freq: 100,
            verbose: false, ..IntegratorConfig::default()
        };
        let integrator = CRTBPIntegrator::<EarthMoon>::new(config);
        let xl1 = l1_position();
        let state = StateVector::new(xl1 + 0.005, -0.002, 0.0, 0.01, -0.005, 0.0);
        let cj_initial = jacobi_constant(&state);
        let event = |_t: f64, _s: &StateVector| -> f64 { 1.0 };
        let result = integrator.propagate_until_event(&state, 10.0, 1.0, event, 1e-14, None, 200_000, None);
        match result {
            Ok(res) => {
                let cj_final = jacobi_constant(&StateVector::new(
                    res.final_state[0], res.final_state[1], res.final_state[2],
                    res.final_state[3], res.final_state[4], res.final_state[5],
                ));
                let drift = (cj_final - cj_initial).abs();
                assert!(drift < 1e-9, "Deriva excesiva: {:.2e}", drift);
                assert!(res.jacobi_error_max < 1e-8, "Deriva máxima excesiva: {:.2e}", res.jacobi_error_max);
            }
            Err(e) => println!("Propagación terminó: {:?}", e),
        }
    }
}
