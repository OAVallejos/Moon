//! low_thrust.rs - Maniobras orbitales lunares con RK45 + control tangencial
//! L-RAIL TUG-BIG: 8x NEXT-C, empuje tangencial retrógrado para descenso

use crate::crtbp::{StateVector, crtbp_derivatives};
use crate::constants::{D_CHAR, T_CHAR, MU_NORMALIZED};

#[derive(Clone)]
pub struct LowThrustConfig {
    pub thrust_max: f64,
    pub isp: f64,
    pub initial_mass: f64,
    pub dry_mass: f64,
    pub max_angle_rate: f64,
    pub min_throttle: f64,
    pub max_throttle: f64,
}

impl Default for LowThrustConfig {
    fn default() -> Self {
        LowThrustConfig {
            thrust_max: 0.64,
            isp: 4190.0,
            initial_mass: 1910.5,
            dry_mass: 1610.5,
            max_angle_rate: 0.1,
            min_throttle: 0.05,
            max_throttle: 1.0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ControlVector {
    pub throttle: f64,
    pub alpha: f64,
    pub beta: f64,
}

impl ControlVector {
    pub fn zero() -> Self {
        ControlVector { throttle: 0.0, alpha: 0.0, beta: 0.0 }
    }
}

pub struct RK45Result {
    pub next_state: StateVector,
    pub h_next: f64,
    pub accepted: bool,
}

pub fn rk45_adaptive_step<F>(
    f: &F,
    t: f64,
    state: &StateVector,
    h: f64,
    tolerance: f64,
) -> RK45Result
where
    F: Fn(f64, &StateVector) -> StateVector,
{
    let k1 = f(t, state);
    let mut s2 = StateVector::zeros();
    for i in 0..6 { s2[i] = state[i] + h * (1.0/4.0) * k1[i]; }
    let k2 = f(t + 0.25 * h, &s2);
    let mut s3 = StateVector::zeros();
    for i in 0..6 { s3[i] = state[i] + h * (3.0/32.0 * k1[i] + 9.0/32.0 * k2[i]); }
    let k3 = f(t + 3.0/8.0 * h, &s3);
    let mut s4 = StateVector::zeros();
    for i in 0..6 { s4[i] = state[i] + h * (1932.0/2197.0 * k1[i] - 7200.0/2197.0 * k2[i] + 7296.0/2197.0 * k3[i]); }
    let k4 = f(t + 12.0/13.0 * h, &s4);
    let mut s5 = StateVector::zeros();
    for i in 0..6 { s5[i] = state[i] + h * (439.0/216.0 * k1[i] - 8.0 * k2[i] + 3680.0/513.0 * k3[i] - 845.0/4104.0 * k4[i]); }
    let k5 = f(t + h, &s5);
    let mut s6 = StateVector::zeros();
    for i in 0..6 { s6[i] = state[i] + h * (-8.0/27.0 * k1[i] + 2.0 * k2[i] - 3544.0/2565.0 * k3[i] + 1859.0/4104.0 * k4[i] - 11.0/40.0 * k5[i]); }
    let k6 = f(t + 0.5 * h, &s6);

    let mut y5 = StateVector::zeros();
    for i in 0..6 {
        y5[i] = state[i] + h * (16.0/135.0 * k1[i] + 6656.0/12825.0 * k3[i] + 28561.0/56430.0 * k4[i] - 9.0/50.0 * k5[i] + 2.0/55.0 * k6[i]);
    }

    let mut error = 0.0;
    for i in 0..6 {
        let err_i = h * (1.0/360.0 * k1[i] - 128.0/4275.0 * k3[i] - 2197.0/75240.0 * k4[i] + 1.0/50.0 * k5[i] + 2.0/55.0 * k6[i]);
        error += err_i * err_i;
    }
    error = error.sqrt();

    let safety = 0.9;
    let scale = if error > 1e-16 { safety * (tolerance / error).powf(0.2) } else { 2.0 };
    let scale = scale.clamp(0.1, 4.0);

    if error <= tolerance {
        RK45Result { next_state: y5, h_next: h * scale, accepted: true }
    } else {
        RK45Result { next_state: *state, h_next: h * scale.clamp(0.1, 0.5), accepted: false }
    }
}

pub fn optimize_low_thrust_transfer(
    initial_state: &[f64; 6],
    target_altitude_km: f64,
    time_of_flight_days: f64,
    config: LowThrustConfig,
) -> Result<(Vec<[f64; 6]>, Vec<ControlVector>), String> {
    let mu = MU_NORMALIZED;
    let moon_x = 1.0 - mu;
    let r_moon_can = 1737.0 * 1000.0 / D_CHAR;
    let tof_canonical = time_of_flight_days * 86400.0 / T_CHAR;
    
    let n_segments = 2000;
    let dt_segment = tof_canonical / n_segments as f64;
    let tolerance = 1e-10;

    let mut state = StateVector::new(
        initial_state[0], initial_state[1], initial_state[2],
        initial_state[3], initial_state[4], initial_state[5],
    );

    let dx0 = state[0] - moon_x;
    let dy0 = state[1];
    let dz0 = state[2];
    let r0 = (dx0*dx0 + dy0*dy0 + dz0*dz0).sqrt();
    let alt0_km = (r0 - r_moon_can) * D_CHAR / 1000.0;

    println!("   Altitud inicial: {:.1} km", alt0_km);
    println!("   Altitud objetivo: {:.1} km", target_altitude_km);

    if alt0_km < 50.0 || alt0_km > 5000.0 {
        return Err(format!("No está en órbita lunar válida (altitud={:.0} km).", alt0_km));
    }

    let climbing = target_altitude_km > alt0_km;
    let delta_alt = (target_altitude_km - alt0_km).abs();
    println!("   Maniobra: {} (Δalt={:.0} km)", if climbing { "ASCENSO" } else { "DESCENSO" }, delta_alt);

    let m_char = 6.045e24;
    let thrust_canonical = config.thrust_max * T_CHAR.powi(2) / (m_char * D_CHAR);
    let ve_m_s = config.isp * 9.80665;
    let m_dot_max = config.thrust_max / ve_m_s;
    let m_dot_canonical = m_dot_max * T_CHAR / m_char;

    let mut mass_canonical = config.initial_mass / m_char;
    let dry_mass_canonical = config.dry_mass / m_char;

    let mut controls = vec![ControlVector::zero(); n_segments];
    let mut trajectory = Vec::with_capacity(n_segments + 1);
    trajectory.push([state[0], state[1], state[2], state[3], state[4], state[5]]);

    for i in 0..n_segments {
        let t_current = i as f64 * dt_segment;

        let dx_moon = state[0] - moon_x;
        let dy_moon = state[1];
        let dz_moon = state[2];
        let r_moon = (dx_moon*dx_moon + dy_moon*dy_moon + dz_moon*dz_moon).sqrt();
        let alt_km = (r_moon - r_moon_can) * D_CHAR / 1000.0;

        let alt_error = target_altitude_km - alt_km;
        let (throttle, alpha, beta) = if alt_error.abs() > 2.0 {
            let throttle = 1.0;
            let vx_rel = state[3];
            let vy_rel = state[4];
            let vz_rel = state[5];
            let v_rel_mag = (vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel).sqrt();

            if v_rel_mag < 1e-10 {
                (0.0, 0.0, 0.0)
            } else {
                let uv_x = vx_rel / v_rel_mag;
                let uv_y = vy_rel / v_rel_mag;
                let uv_z = vz_rel / v_rel_mag;
                let sign = if climbing { 1.0 } else { -1.0 };
                
                let dir_x = sign * uv_x;
                let dir_y = sign * uv_y;
                let dir_z = sign * uv_z;

                let alpha = dir_y.atan2(dir_x);
                let beta = dir_z.clamp(-1.0, 1.0).asin();
                (throttle, alpha, beta)
            }
        } else {
            (0.0, 0.0, 0.0)
        };

        controls[i] = ControlVector { throttle, alpha, beta };

        let thrust_dir_x = alpha.cos() * beta.cos();
        let thrust_dir_y = alpha.sin() * beta.cos();
        let thrust_dir_z = beta.sin();

        let thrust_acc = thrust_canonical / mass_canonical;
        let acc_thrust_x = throttle * thrust_acc * thrust_dir_x;
        let acc_thrust_y = throttle * thrust_acc * thrust_dir_y;
        let acc_thrust_z = throttle * thrust_acc * thrust_dir_z;

        let derivatives = |t: f64, s: &StateVector| -> StateVector {
            let mut deriv = crtbp_derivatives(t, s);
            deriv[3] += acc_thrust_x;
            deriv[4] += acc_thrust_y;
            deriv[5] += acc_thrust_z;
            deriv
        };

        let mut t_sub = t_current;
        let target_t = t_current + dt_segment;
        let mut h_current = dt_segment * 0.1;
        let mut max_substeps = 0;

        while t_sub < target_t && max_substeps < 10000 {
            if t_sub + h_current > target_t { h_current = target_t - t_sub; }
            let res = rk45_adaptive_step(&derivatives, t_sub, &state, h_current, tolerance);
            if res.accepted {
                state = res.next_state;
                t_sub += h_current;
            }
            h_current = res.h_next;
            max_substeps += 1;
        }

        mass_canonical -= throttle * m_dot_canonical * dt_segment;
        mass_canonical = mass_canonical.max(dry_mass_canonical);
        trajectory.push([state[0], state[1], state[2], state[3], state[4], state[5]]);
    }

    let final_state = trajectory.last().unwrap();
    let dx = final_state[0] - moon_x;
    let dy = final_state[1];
    let dz = final_state[2];
    let r_final = (dx*dx + dy*dy + dz*dz).sqrt();
    let alt_final_km = (r_final - r_moon_can) * D_CHAR / 1000.0;
    let xenon_used_kg = (config.initial_mass - mass_canonical * m_char).max(0.0);

    println!("   Altitud final: {:.1} km (objetivo: {:.1} km)", alt_final_km, target_altitude_km);
    println!("   Xenón consumido: {:.1} kg", xenon_used_kg);

    let tolerance_km = 50.0;
    let alt_error_final = (alt_final_km - target_altitude_km).abs();

    if alt_error_final > tolerance_km {
        return Err(format!("No se alcanzó la órbita objetivo. Altitud final: {:.1} km, objetivo: {:.1} km, error: {:.1} km", alt_final_km, target_altitude_km, alt_error_final));
    }

    println!("   ✅ Maniobra exitosa (error: {:.1} km)", alt_error_final);
    Ok((trajectory, controls))
}
