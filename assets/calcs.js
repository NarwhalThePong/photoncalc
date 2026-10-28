/* PhysCalc calculator math. Pure functions, no DOM. Works in browser (window.PC) and Node (module.exports). */
(function (root) {
  'use strict';
  var C = 299792458;            // m/s, exact (SI)
  var H = 6.62607015e-34;       // J s, exact (SI 2019)
  var QE = 1.602176634e-19;     // C, exact (SI 2019)
  var PI = Math.PI;
  var D2R = PI / 180, R2D = 180 / PI;
  function req(cond, msg) { if (!cond) throw new Error(msg); }
  function pos(x, name) { req(isFinite(x) && x > 0, name + ' must be a positive number'); }

  // 1. wavelength / frequency / photon energy. unit: nm, um, THz, eV, cm-1
  function photon(v, unit) {
    pos(v, 'Value');
    var lam; // metres
    if (unit === 'nm') lam = v * 1e-9;
    else if (unit === 'um') lam = v * 1e-6;
    else if (unit === 'THz') lam = C / (v * 1e12);
    else if (unit === 'eV') lam = H * C / (v * QE);
    else if (unit === 'cm-1') lam = 1 / (v * 100);
    else throw new Error('Unknown unit');
    var f = C / lam, EJ = H * f;
    return { lambda_nm: lam * 1e9, lambda_um: lam * 1e6, f_THz: f / 1e12, E_eV: EJ / QE, E_J: EJ, wavenumber_cm: 1 / (lam * 100) };
  }

  // 2. optical / RF power: unit dBm, dBW, mW, W
  function power(v, unit) {
    var mW;
    if (unit === 'dBm') { req(isFinite(v), 'Value must be a number'); mW = Math.pow(10, v / 10); }
    else if (unit === 'dBW') { req(isFinite(v), 'Value must be a number'); mW = 1000 * Math.pow(10, v / 10); }
    else if (unit === 'mW') { pos(v, 'Power'); mW = v; }
    else if (unit === 'W') { pos(v, 'Power'); mW = v * 1000; }
    else throw new Error('Unknown unit');
    var dBm = 10 * Math.log10(mW);
    return { mW: mW, W: mW / 1000, uW: mW * 1000, dBm: dBm, dBW: dBm - 30 };
  }
  function dbRatio(db) { req(isFinite(db), 'dB must be a number'); return { power_ratio: Math.pow(10, db / 10), amplitude_ratio: Math.pow(10, db / 20) }; }

  // 3. Gaussian beam. lambda nm, w0 um (1/e^2 radius), M2, z mm
  function gaussian(lam_nm, w0_um, M2, z_mm) {
    pos(lam_nm, 'Wavelength'); pos(w0_um, 'Waist radius'); req(isFinite(M2) && M2 >= 1, 'M² must be ≥ 1'); req(isFinite(z_mm), 'Distance must be a number');
    var lam = lam_nm * 1e-9, w0 = w0_um * 1e-6, z = z_mm * 1e-3;
    var zR = PI * w0 * w0 / (M2 * lam);
    var theta = M2 * lam / (PI * w0);
    var wz = w0 * Math.sqrt(1 + (z / zR) * (z / zR));
    var Rz = z === 0 ? Infinity : z * (1 + (zR / z) * (zR / z));
    return { zR_mm: zR * 1e3, theta_half_mrad: theta * 1e3, theta_full_mrad: 2 * theta * 1e3, w_z_um: wz * 1e6, d_z_um: 2 * wz * 1e6, R_z_mm: Rz * 1e3, BPP_mm_mrad: w0 * 1e3 * theta * 1e3, confocal_mm: 2 * zR * 1e3, gouy_deg: Math.atan(z / zR) * R2D };
  }

  // 4. thin lens: f, object distance (same units). Real-is-positive convention 1/f = 1/do + 1/di
  function thinLens(f, dO) {
    req(isFinite(f) && f !== 0, 'Focal length must be non-zero'); pos(dO, 'Object distance');
    if (Math.abs(dO - f) < 1e-12 * Math.abs(f)) return { di: Infinity, M: Infinity, note: 'Object at the focal point: image at infinity (collimated output).' };
    var di = 1 / (1 / f - 1 / dO);
    var M = -di / dO;
    return { di: di, M: M, image: di > 0 ? 'real, inverted' : 'virtual, upright', power_D_if_mm: 1000 / f };
  }

  // 5. numerical aperture from indices
  function numAperture(n1, n2, n0) {
    pos(n1, 'Core index'); pos(n2, 'Cladding index'); pos(n0, 'Outside index');
    req(n1 > n2, 'Core index must be larger than cladding index');
    var NA = Math.sqrt(n1 * n1 - n2 * n2);
    var s = NA / n0;
    var acc = s >= 1 ? 90 : Math.asin(s) * R2D;
    return { NA: NA, accept_half_deg: acc, accept_full_deg: 2 * acc, delta: (n1 - n2) / n1, critical_deg: Math.asin(n2 / n1) * R2D };
  }

  // 6. fiber V-number. a core radius um, NA, lambda nm
  var J0ZERO = 2.404825557695773;
  function vNumber(a_um, NA, lam_nm) {
    pos(a_um, 'Core radius'); pos(NA, 'NA'); pos(lam_nm, 'Wavelength');
    var V = 2 * PI * a_um * 1e3 * NA / lam_nm;
    // Marcuse mode-field radius estimate, valid ~0.8 < V < 2.5
    var w = a_um * (0.65 + 1.619 / Math.pow(V, 1.5) + 2.879 / Math.pow(V, 6));
    return { V: V, single_mode: V < J0ZERO, cutoff_nm: 2 * PI * a_um * 1e3 * NA / J0ZERO, modes_approx: V * V / 2, mfd_um: 2 * w };
  }

  // 7. Snell / critical / Brewster / Fresnel
  function snell(n1, n2, thi_deg) {
    pos(n1, 'n1'); pos(n2, 'n2'); req(isFinite(thi_deg) && thi_deg >= 0 && thi_deg < 90, 'Incidence angle must be 0–89.999°');
    var ti = thi_deg * D2R, s = n1 * Math.sin(ti) / n2;
    var out = { brewster_deg: Math.atan(n2 / n1) * R2D, critical_deg: n1 > n2 ? Math.asin(n2 / n1) * R2D : null };
    if (s > 1) { out.tir = true; out.theta_t_deg = null; out.Rs = 1; out.Rp = 1; return out; }
    var tt = Math.asin(s), ci = Math.cos(ti), ct = Math.cos(tt);
    var rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
    var rp = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct);
    out.tir = false; out.theta_t_deg = tt * R2D; out.Rs = rs * rs; out.Rp = rp * rp; out.Ravg = (rs * rs + rp * rp) / 2;
    return out;
  }

  // 8. pulsed laser. Pavg W, frep Hz, tau FWHM s, shape, beam 1/e^2 diameter um, lambda nm
  var SHAPE_K = { gaussian: 2 * Math.sqrt(Math.LN2 / PI), sech2: Math.log(1 + Math.SQRT2), square: 1 };
  function pulse(Pavg, frep, tau, shape, D_um, lam_nm) {
    pos(Pavg, 'Average power'); pos(frep, 'Repetition rate'); pos(tau, 'Pulse duration'); pos(D_um, 'Beam diameter'); pos(lam_nm, 'Wavelength');
    var k = SHAPE_K[shape]; req(k, 'Unknown pulse shape');
    var E = Pavg / frep, Pk = k * E / tau, w_cm = D_um * 1e-4 / 2;
    return { E_J: E, Ppeak_W: Pk, duty: tau * frep, fluence_Jcm2: 2 * E / (PI * w_cm * w_cm), intensity_Wcm2: 2 * Pk / (PI * w_cm * w_cm), photons: E * lam_nm * 1e-9 / (H * C), k: k };
  }

  // 9. PPLN QPM poling period, 5 mol% MgO:congruent LN, extraordinary index (Gayer et al., Appl. Phys. B 91, 343 (2008); erratum 101, 481 (2010))
  var G = { a: [5.756, 0.0983, 0.2020, 189.32, 12.52, 1.32e-2], b: [2.860e-6, 4.700e-8, 6.113e-8, 1.516e-4] };
  function neMgLN(lam_um, T) {
    var a = G.a, b = G.b, f = (T - 24.5) * (T + 570.82), l2 = lam_um * lam_um, a3 = a[2] + b[2] * f;
    return Math.sqrt(a[0] + b[0] * f + (a[1] + b[1] * f) / (l2 - a3 * a3) + (a[3] + b[3] * f) / (l2 - a[4] * a[4]) - a[5] * l2);
  }
  function ppln(lp_nm, T, m) {
    pos(lp_nm, 'Pump wavelength'); req(isFinite(T) && T >= 0 && T <= 250, 'Temperature should be 0–250 °C'); req(Number.isInteger(m) && m >= 1, 'QPM order must be a positive integer');
    req(lp_nm >= 700 && lp_nm <= 5000, 'Pump wavelength outside 700–5000 nm (SH must stay inside the fit range)');
    var lp = lp_nm / 1000, ls = lp / 2, np = neMgLN(lp, T), ns = neMgLN(ls, T);
    var L1 = 1 / (ns / ls - 2 * np / lp); // um
    return { n_pump: np, n_sh: ns, lambda_sh_nm: lp_nm / 2, period_um: m * L1, period_first_um: L1, coherence_len_um: L1 / 2 };
  }

  // 10. Fabry-Perot. L mm, n, R1, R2, lambda nm
  function fabryPerot(L_mm, n, R1, R2, lam_nm) {
    pos(L_mm, 'Length'); pos(n, 'Index'); pos(lam_nm, 'Wavelength');
    req(R1 > 0 && R1 < 1 && R2 > 0 && R2 < 1, 'Reflectivities must be between 0 and 1');
    var L = L_mm * 1e-3, fsr = C / (2 * n * L), Rg = Math.sqrt(R1 * R2);
    var F = PI * Math.pow(R1 * R2, 0.25) / (1 - Rg);
    var lam = lam_nm * 1e-9, nu = C / lam;
    return { fsr_Hz: fsr, fsr_nm: lam * lam / (2 * n * L) * 1e9, finesse: F, fwhm_Hz: fsr / F, Q: nu / (fsr / F), coef_F: 4 * Rg / Math.pow(1 - Rg, 2) };
  }

  // 11. grating equation  m*lambda = d (sin a + sin b), angles from normal, same-side positive
  function grating(gpmm, lam_nm, alpha_deg, m) {
    pos(gpmm, 'Groove density'); pos(lam_nm, 'Wavelength'); req(isFinite(alpha_deg) && Math.abs(alpha_deg) < 90, 'Incidence angle must be between −90 and 90°'); req(Number.isInteger(m), 'Order must be an integer');
    var d_nm = 1e6 / gpmm, s = m * lam_nm / d_nm - Math.sin(alpha_deg * D2R);
    var lit = lam_nm / (2 * d_nm);
    var out = { d_um: d_nm / 1000, littrow_deg: lit <= 1 ? Math.asin(lit) * R2D : null, max_order: Math.floor(2 * d_nm / lam_nm) };
    if (Math.abs(s) > 1) { out.evanescent = true; out.beta_deg = null; return out; }
    var b = Math.asin(s); out.evanescent = false; out.beta_deg = b * R2D;
    out.dispersion_deg_per_nm = (m / (d_nm * Math.cos(b))) * R2D; return out;
  }

  // 12. LC resonance (+ optional series R for Q, + RC cutoff)
  function lc(L_H, C_F, R_ohm) {
    pos(L_H, 'Inductance'); pos(C_F, 'Capacitance');
    var w0 = 1 / Math.sqrt(L_H * C_F), out = { f0_Hz: w0 / (2 * PI), w0: w0, Z0: Math.sqrt(L_H / C_F) };
    if (isFinite(R_ohm) && R_ohm > 0) { out.Q_series = Math.sqrt(L_H / C_F) / R_ohm; out.bw_Hz = out.f0_Hz / out.Q_series; out.rc_cutoff_Hz = 1 / (2 * PI * R_ohm * C_F); }
    return out;
  }

  // 13. general three-wave QPM in MgO:LN (SFG/DFG/OPO), pump > signal >= idler frequency. lp, ls nm
  function qpm3(lp_nm, ls_nm, T, m) {
    pos(lp_nm, 'Pump wavelength'); pos(ls_nm, 'Signal wavelength'); req(isFinite(T) && T >= 0 && T <= 250, 'Temperature should be 0–250 °C'); req(Number.isInteger(m) && m >= 1, 'QPM order must be a positive integer');
    req(ls_nm > lp_nm, 'Signal wavelength must be longer than the pump');
    var li_nm = 1 / (1 / lp_nm - 1 / ls_nm);
    req(lp_nm >= 400 && li_nm <= 5000, 'All three wavelengths must stay within 400–5000 nm (Sellmeier fit range)');
    var lp = lp_nm / 1000, ls = ls_nm / 1000, li = li_nm / 1000;
    var np = neMgLN(lp, T), ns = neMgLN(ls, T), ni = neMgLN(li, T);
    var L1 = 1 / (np / lp - ns / ls - ni / li);
    return { lambda_idler_nm: li_nm, n_pump: np, n_signal: ns, n_idler: ni, period_um: m * L1, period_first_um: L1 };
  }
  // 14. SHG phase-matching temperature for a given poling period (bisection over 0–250 C)
  function pplnTemp(lp_nm, period_um, m) {
    pos(period_um, 'Poling period');
    var f = function (T) { return ppln(lp_nm, T, m).period_um - period_um; };
    var lo = 0, hi = 250, flo = f(lo), fhi = f(hi);
    if (flo * fhi > 0) return { T_C: null, period_at_0: flo + period_um, period_at_250: fhi + period_um };
    for (var i = 0; i < 80; i++) { var mid = (lo + hi) / 2, fm = f(mid); if (fm * flo <= 0) { hi = mid; } else { lo = mid; flo = fm; } }
    var T = (lo + hi) / 2, d = (ppln(lp_nm, T + 0.5, m).period_um - ppln(lp_nm, T - 0.5, m).period_um);
    return { T_C: T, dLambda_dT_nm_per_K: d * 1000, period_at_0: f(0) + period_um, period_at_250: f(250) + period_um };
  }

  // 15. microring resonator. R um, ng, neff, lambda nm, loss dB/cm, power coupling kappa^2 (each coupler), type 'allpass'|'adddrop'
  function ring(R_um, ng, lam_nm, loss_dBcm, k2, type) {
    pos(R_um, 'Radius'); pos(ng, 'Group index'); pos(lam_nm, 'Wavelength'); req(isFinite(loss_dBcm) && loss_dBcm >= 0, 'Loss must be ≥ 0 dB/cm');
    req(k2 > 0 && k2 < 1, 'Power coupling κ² must be between 0 and 1'); req(type === 'allpass' || type === 'adddrop', 'Unknown configuration');
    var L = 2 * PI * R_um * 1e-6, lam = lam_nm * 1e-9;
    var alpha = loss_dBcm * 100 / (10 / Math.LN10);   // power attenuation, 1/m
    var a = Math.exp(-alpha * L / 2), r = Math.sqrt(1 - k2);
    var ra = type === 'allpass' ? r * a : r * r * a;
    var fsr_m = lam * lam / (ng * L), fwhm_m = (1 - ra) * lam * lam / (PI * ng * L * Math.sqrt(ra));
    var out = { L_um: L * 1e6, fsr_nm: fsr_m * 1e9, fsr_GHz: C / (ng * L) / 1e9, fwhm_pm: fwhm_m * 1e12, finesse: fsr_m / fwhm_m, Q_loaded: lam / fwhm_m,
      Q_intrinsic: alpha > 0 ? 2 * PI * ng / (lam * alpha) : Infinity, roundtrip_loss_dB: loss_dBcm * L * 100, a: a, r: r };
    if (type === 'allpass') { var Tm = (a - r) * (a - r) / ((1 - r * a) * (1 - r * a)); out.T_min = Tm; out.ER_dB = Tm > 0 ? -10 * Math.log10(Tm) : Infinity; out.regime = Math.abs(r - a) < 1e-9 ? 'critically coupled' : (r < a ? 'over-coupled' : 'under-coupled'); }
    else { out.T_through_min = Math.pow(r - r * a, 2) / Math.pow(1 - r * r * a, 2); out.T_drop_max = Math.pow(1 - r * r, 2) * a / Math.pow(1 - r * r * a, 2); }
    return out;
  }

  // 16. Sellmeier index, group index, GVD. material key, lambda nm
  var SELL = {
    fused_silica: { B: [0.6961663, 0.4079426, 0.8974794], Cs: [0.0684043, 0.1162414, 9.896161], rng: [0.21, 3.71] },   // Malitson 1965
    nbk7: { B: [1.03961212, 0.231792344, 1.01046945], C: [0.00600069867, 0.0200179144, 103.560653], rng: [0.3, 2.5] },    // SCHOTT
    sapphire_o: { B: [1.4313493, 0.65054713, 5.3414021], Cs: [0.0726631, 0.1193242, 18.028251], rng: [0.2, 5.0] },       // Malitson & Dodge 1972
    silicon: { fn: function (l) { var l1 = 1.1071, l2 = l * l; return Math.sqrt(11.6858 + 0.939816 / l2 + 0.00810461 * l1 * l1 / (l2 - l1 * l1)); }, rng: [1.36, 11] }, // Salzberg & Villa 1957
    mgln_e: { fn: function (l) { return neMgLN(l, 24.5); }, rng: [0.5, 4.0] }                                             // Gayer 2008, 24.5 C
  };
  function sellN(mat, l) {
    var M = SELL[mat]; if (M.fn) return M.fn(l);
    var l2 = l * l, n2 = 1; for (var i = 0; i < 3; i++) { var c = M.C ? M.C[i] : M.Cs[i] * M.Cs[i]; n2 += M.B[i] * l2 / (l2 - c); } return Math.sqrt(n2);
  }
  function sellmeier(mat, lam_nm) {
    var M = SELL[mat]; req(M, 'Unknown material'); pos(lam_nm, 'Wavelength');
    var l = lam_nm / 1000; req(l >= M.rng[0] && l <= M.rng[1], 'Wavelength outside this fit\'s range (' + M.rng[0] * 1000 + '–' + M.rng[1] * 1000 + ' nm)');
    var h = 1e-3 * l, n = sellN(mat, l), np = sellN(mat, l + h), nm = sellN(mat, l - h);
    var d1 = (np - nm) / (2 * h), d2 = (np - 2 * n + nm) / (h * h); // per um, per um^2
    var ng = n - l * d1;
    var b2 = Math.pow(l * 1e-6, 3) / (2 * PI * C * C) * d2 * 1e12; // s^2/m
    return { n: n, ng: ng, dn_dlam_per_um: d1, gvd_fs2_per_mm: b2 * 1e30 / 1e3, D_ps_nm_km: -(l * 1e-6) / C * d2 * 1e12 * 1e12 / 1e9 * 1e3, range_nm: [M.rng[0] * 1000, M.rng[1] * 1000] };
  }

  // 17. Gaussian mode coupling efficiency: waists w1,w2 um, lateral offset d um, tilt deg, lambda nm, gap medium index n
  function coupling(w1, w2, d_um, tilt_deg, lam_nm, n) {
    pos(w1, 'Mode radius 1'); pos(w2, 'Mode radius 2'); pos(lam_nm, 'Wavelength'); pos(n, 'Index');
    req(isFinite(d_um) && d_um >= 0, 'Offset must be ≥ 0'); req(isFinite(tilt_deg) && tilt_deg >= 0 && tilt_deg < 90, 'Tilt must be 0–89°');
    var s = w1 * w1 + w2 * w2, lam = lam_nm / 1000, th = tilt_deg * D2R;
    var eS = Math.pow(2 * w1 * w2 / s, 2), eD = Math.exp(-2 * d_um * d_um / s), eT = Math.exp(-2 * PI * PI * n * n * th * th * w1 * w1 * w2 * w2 / (lam * lam * s));
    var eta = eS * eD * eT;
    var dB = function (x) { return -10 * Math.log10(x) + 0; };
    return { eta: eta, loss_dB: dB(eta), size_dB: dB(eS), offset_dB: dB(eD), tilt_dB: dB(eT) };
  }

  // 18. grating coupler phase matching: lambda nm, n_eff, period nm, cladding index, order m, target angle deg (in cladding)
  function gratingCoupler(lam_nm, neff, per_nm, nc, m, target_deg) {
    pos(lam_nm, 'Wavelength'); pos(neff, 'Effective index'); pos(per_nm, 'Grating period'); pos(nc, 'Cladding index'); req(m >= 1, 'Order must be ≥ 1');
    req(isFinite(target_deg) && Math.abs(target_deg) < 90, 'Target angle must be between −90° and 90°');
    var s = (neff - m * lam_nm / per_nm) / nc, sa = s * nc;   // sin(theta) in cladding; n_c sin = air-side sin (Snell into air)
    var out = { sin_clad: s, theta_clad_deg: Math.abs(s) <= 1 ? Math.asin(s) / D2R : null, theta_air_deg: Math.abs(sa) <= 1 ? Math.asin(sa) / D2R : null,
      period_for_target_nm: m * lam_nm / (neff - nc * Math.sin(target_deg * D2R)), dtheta_dlam_deg_per_nm: null };
    if (out.theta_clad_deg !== null) out.dtheta_dlam_deg_per_nm = -m / (per_nm * nc * Math.cos(out.theta_clad_deg * D2R)) / D2R; // fixed n_eff (no waveguide dispersion)
    return out;
  }

  // 19. time-bandwidth product. tau fs (FWHM intensity), bandwidth nm (FWHM), centre lambda nm, shape 'gauss'|'sech2', GDD fs^2
  var TBP = { gauss: 2 * Math.LN2 / PI, sech2: 0.314833 };
  function tbp(tau_fs, dl_nm, lam_nm, shape, gdd_fs2) {
    pos(tau_fs, 'Pulse duration'); pos(dl_nm, 'Bandwidth'); pos(lam_nm, 'Wavelength'); req(TBP[shape], 'Unknown pulse shape'); req(isFinite(gdd_fs2), 'GDD must be a number');
    var k = TBP[shape], lam = lam_nm * 1e-9, dnu = C * dl_nm * 1e-9 / (lam * lam), tau = tau_fs * 1e-15;
    var tl = k / dnu, dnu_tl = k / tau, out = { k: k, dnu_THz: dnu / 1e12, tbp: dnu * tau, tl_fs: tl * 1e15, ratio: tau / tl,
      dl_tl_nm: dnu_tl * lam * lam / C * 1e9, dnu_tl_THz: dnu_tl / 1e12, gdd_out_fs: null };
    if (shape === 'gauss') { var t0 = tl * 1e15; out.gdd_out_fs = t0 * Math.sqrt(1 + Math.pow(4 * Math.LN2 * gdd_fs2 / (t0 * t0), 2)); } // TL Gaussian at this bandwidth after GDD
    return out;
  }

  // 20. photodiode: lambda nm, external QE (0-1), optical power W, load ohm
  function photodiode(lam_nm, qe, P_W, RL) {
    pos(lam_nm, 'Wavelength'); req(qe > 0 && qe <= 1, 'Quantum efficiency must be between 0 and 100 %'); req(isFinite(P_W) && P_W >= 0, 'Power must be ≥ 0'); req(isFinite(RL) && RL > 0, 'Load must be > 0 Ω');
    var Rmax = QE * lam_nm * 1e-9 / (H * C), R = qe * Rmax, I = R * P_W;
    return { R_AW: R, Rmax_AW: Rmax, I_A: I, V_V: I * RL, shot_A_rtHz: Math.sqrt(2 * QE * I), photon_rate: P_W * lam_nm * 1e-9 / (H * C) };
  }

  // 21. focused spot / beam expander. lambda nm, input 1/e^2 diameter mm, expander magnification, focal length mm, M^2
  function focusSpot(lam_nm, D_mm, mag, f_mm, m2) {
    pos(lam_nm, 'Wavelength'); pos(D_mm, 'Beam diameter'); pos(mag, 'Expander magnification'); pos(f_mm, 'Focal length'); req(isFinite(m2) && m2 >= 1, 'M² must be ≥ 1');
    var lam = lam_nm * 1e-6, De = D_mm * mag, w = De / 2, w0 = m2 * lam * f_mm / (PI * w), zR = PI * w0 * w0 / (m2 * lam);
    return { D_exp_mm: De, w0_um: w0 * 1e3, spot_um: 2 * w0 * 1e3, zR_um: zR * 1e3, dof_um: 2 * zR * 1e3, NA: w / f_mm, fnum: f_mm / De,
      theta_coll_urad: m2 * lam / (PI * w) * 1e6, paraxial_ok: w / f_mm <= 0.25 };
  }

  // 22. SHG coherence length. material key (SELL), pump lambda nm, crystal length mm (0 = skip)
  function cohLen(mat, lp_nm, L_mm) {
    var M = SELL[mat]; req(M, 'Unknown material'); pos(lp_nm, 'Pump wavelength'); req(isFinite(L_mm) && L_mm >= 0, 'Crystal length must be ≥ 0');
    var lp = lp_nm / 1000, ls = lp / 2;
    req(ls >= M.rng[0] && lp <= M.rng[1], 'Pump and second harmonic must both lie inside this fit\'s range (' + M.rng[0] * 1000 + '–' + M.rng[1] * 1000 + ' nm), so the pump must be ' + M.rng[0] * 2000 + '–' + M.rng[1] * 1000 + ' nm');
    var n1 = sellN(mat, lp), n2 = sellN(mat, ls), dn = n2 - n1;
    req(Math.abs(dn) > 1e-12, 'Indices are equal: already phase matched (infinite coherence length)');
    var dk = 4 * PI * dn / lp; // 1/um, dk = k(2w) - 2k(w)
    var Lc = PI / Math.abs(dk); // um, = lp/(4|dn|)
    var out = { n_pump: n1, n_sh: n2, dn: dn, dk_per_mm: dk * 1e3, Lc_um: Lc, period1_um: 2 * Lc, lambda_sh_nm: lp_nm / 2, sinc2: null, n_coh: null };
    if (L_mm > 0) { var x = dk * L_mm * 1e3 / 2; out.sinc2 = Math.pow(Math.sin(x) / x, 2); out.n_coh = L_mm * 1e3 / Lc; }
    return out;
  }


  // 23. TFLN ridge-waveguide QPM poling period (type-0 SHG, TE00 in x-cut / extraordinary). Effective index method (EIM):
  // vertical asymmetric slab (TE) for the ridge and the etched sides, then a lateral symmetric slab (TM) of width w.
  // lp nm, film thickness t nm, etch depth e nm, top width w nm, material 'cln'|'mgln', cover 'air'|'sio2', T C, order m
  function neCLN(lam_um, T) { // congruent LN, extraordinary (Jundt, Opt. Lett. 22, 1553 (1997))
    var f = (T - 24.5) * (T + 570.82), l2 = lam_um * lam_um, a3 = 0.20692 - 0.89e-8 * f;
    return Math.sqrt(5.35583 + 4.629e-7 * f + (0.100473 + 3.862e-8 * f) / (l2 - a3 * a3) + (100 + 2.657e-5 * f) / (l2 - 11.34927 * 11.34927) - 1.5334e-2 * l2);
  }
  // fundamental mode of a 3-layer slab: core nf, thickness d (um), substrate ns, cover nc; tm = TM equation. Returns N or null if cut off.
  function slabN(k0, d, nf, ns, nc, tm) {
    var lo = Math.max(ns, nc), F = function (N) {
      var kf = Math.sqrt(nf * nf - N * N), gs = Math.sqrt(N * N - ns * ns), gc = Math.sqrt(N * N - nc * nc);
      var rs = tm ? nf * nf / (ns * ns) : 1, rc = tm ? nf * nf / (nc * nc) : 1;
      return k0 * d * kf - Math.atan(rs * gs / kf) - Math.atan(rc * gc / kf);
    };
    var a = lo + 1e-12, b = nf - 1e-12;
    if (F(a) <= 0) return null; // below cutoff
    for (var i = 0; i < 200; i++) { var m = (a + b) / 2; if (F(m) > 0) a = m; else b = m; }
    return (a + b) / 2;
  }
  function tflnNeff(lam_um, t, e, w, nf, ns, nc) {
    var k0 = 2 * PI / lam_um, Nr = slabN(k0, t, nf, ns, nc, false);
    req(Nr !== null, 'The film is too thin to guide at ' + (lam_um * 1000).toFixed(0) + ' nm');
    var Nside = e >= t ? nc : slabN(k0, t - e, nf, ns, nc, false); if (Nside === null) Nside = ns > nc ? ns : nc;
    var N = slabN(k0, w, Nr, Nside, Nside, true);
    return { ridge: Nr, side: Nside, neff: N === null ? Nside : N };
  }
  function tfln(lp_nm, t_nm, e_nm, w_nm, mat, cover, T, m) {
    pos(lp_nm, 'Pump wavelength'); pos(t_nm, 'Film thickness'); pos(w_nm, 'Ridge width'); req(isFinite(e_nm) && e_nm >= 0 && e_nm <= t_nm, 'Etch depth must be between 0 and the film thickness');
    req(isFinite(T) && T >= 0 && T <= 200, 'Temperature should be 0–200 °C'); req(Number.isInteger(m) && m >= 1, 'QPM order must be a positive integer');
    req(lp_nm >= 1000 && lp_nm <= 3000, 'Pump wavelength outside 1000–3000 nm (both waves must stay in the LN and SiO₂ fit ranges)');
    req(mat === 'cln' || mat === 'mgln', 'Unknown film material'); req(cover === 'air' || cover === 'sio2', 'Unknown cladding');
    var lp = lp_nm / 1000, ls = lp / 2, t = t_nm / 1000, e = e_nm / 1000, w = w_nm / 1000;
    var nL = function (l) { return mat === 'cln' ? neCLN(l, T) : neMgLN(l, T); };
    var nS = function (l) { return sellN('fused_silica', l); };
    var P = tflnNeff(lp, t, e, w, nL(lp), nS(lp), cover === 'air' ? 1 : nS(lp));
    var S = tflnNeff(ls, t, e, w, nL(ls), nS(ls), cover === 'air' ? 1 : nS(ls));
    req(S.neff > P.neff, 'No positive mismatch: check the geometry');
    var L1 = lp / (2 * (S.neff - P.neff)), Lb = 1 / (nL(ls) / ls - 2 * nL(lp) / lp);
    return { neff_pump: P.neff, neff_sh: S.neff, n_film_pump: nL(lp), n_film_sh: nL(ls), ridge_pump: P.ridge, side_pump: P.side,
      period_um: m * L1, period_first_um: L1, period_bulk_um: Lb, lambda_sh_nm: lp_nm / 2, lateral_guided: P.neff > P.side };
  }

  // 24. electro-optic (Pockels) modulator half-wave voltage. lambda nm, electrode gap g um, electrode length L mm,
  // material 'tfln'|'bulkln' (congruent LN n_e from Jundt at 25 C, r33 = 30.9 pm/V) or 'custom' (n, r pm/V), overlap Gamma (0-1],
  // config 'pp' (push-pull MZM, Vpi halved) | 'single' (phase modulator or single-arm MZM). Vpi = lambda g / (n^3 r Gamma L)
  var R33LN = 30.9; // pm/V, LN r33 at 1550 nm (low-frequency/unclamped value commonly used in TFLN modulator papers)
  function eoVpi(lam_nm, g_um, L_mm, mat, n, r_pmV, gamma, config) {
    pos(lam_nm, 'Wavelength'); pos(g_um, 'Electrode gap'); pos(L_mm, 'Electrode length');
    req(isFinite(gamma) && gamma > 0 && gamma <= 1, 'Overlap factor Γ must be between 0 and 1');
    req(config === 'pp' || config === 'single', 'Unknown modulator configuration');
    var nn, rr;
    if (mat === 'tfln' || mat === 'bulkln') { req(lam_nm >= 400 && lam_nm <= 5000, 'Wavelength outside 400–5000 nm (LN Sellmeier fit range)'); nn = neCLN(lam_nm / 1000, 25); rr = R33LN; }
    else if (mat === 'custom') { pos(n, 'Refractive index n'); pos(r_pmV, 'Electro-optic coefficient r'); nn = n; rr = r_pmV; }
    else throw new Error('Unknown material');
    var n3r = nn * nn * nn * rr * 1e-12;                                         // m/V
    var vpm = lam_nm * 1e-9 * g_um * 1e-6 / (n3r * gamma * L_mm * 1e-3);          // V, single-arm phase modulator
    var vpi = config === 'pp' ? vpm / 2 : vpm;
    return { n: nn, r_pmV: rr, n3r_pmV: n3r * 1e12, vpi_pm_V: vpm, vpi_V: vpi, vpiL_Vcm: vpi * L_mm / 10, vpmL_Vcm: vpm * L_mm / 10,
      dphi_rad_per_V: PI / vpi, dphi_deg_per_V: 180 / vpi, push_pull: config === 'pp' };
  }

  // 25. laser coherence length from linewidth. lambda nm, linewidth value + unit (Hz/kHz/MHz/GHz/pm/nm), lineshape, medium (group) index n,
  // optional path difference opd_m (medium length, 0 = skip). tau_c: Lorentzian 1/(pi dnu), Gaussian sqrt(2 ln2/pi)/dnu, simple 1/dnu.
  function coherence(lam_nm, w, unit, shape, n, opd_m) {
    pos(lam_nm, 'Wavelength'); pos(w, 'Linewidth'); pos(n, 'Group index n');
    req(isFinite(opd_m) && opd_m >= 0, 'Path difference must be zero or positive');
    var lam = lam_nm * 1e-9, fu = { Hz: 1, kHz: 1e3, MHz: 1e6, GHz: 1e9 }, lu = { pm: 1e-12, nm: 1e-9 }, dnu;
    if (fu[unit]) dnu = w * fu[unit]; else if (lu[unit]) { req(w * lu[unit] < lam, 'Linewidth must be smaller than the wavelength'); dnu = C * w * lu[unit] / (lam * lam); } else throw new Error('Unknown linewidth unit');
    var tau;
    if (shape === 'lorentz') tau = 1 / (PI * dnu); else if (shape === 'gauss') tau = Math.sqrt(2 * Math.LN2 / PI) / dnu; else if (shape === 'simple') tau = 1 / dnu; else throw new Error('Unknown lineshape');
    var vis = null;
    if (opd_m > 0) { var t = n * opd_m / C; vis = shape === 'lorentz' ? Math.exp(-PI * dnu * t) : shape === 'gauss' ? Math.exp(-Math.pow(PI * dnu * t, 2) / (4 * Math.LN2)) : null; }
    return { dnu_Hz: dnu, dlam_pm: dnu * lam * lam / C * 1e12, tau_s: tau, Lc_vac_m: C * tau, Lc_m: C * tau / n, Lc_simple_m: C / (n * dnu), visibility: vis, shape: shape };
  }

  // 26. fiber chromatic dispersion. lambda nm, D ps/(nm km), L km, source spectral width dlam nm (0 = skip),
  // input transform-limited Gaussian pulse FWHM ps (0 = skip). beta2 = -D lambda^2/(2 pi c); T0 = T_fwhm/(2 sqrt(ln2)); L_D = T0^2/|beta2|.
  function dispersion(lam_nm, D, L_km, dlam_nm, tfwhm_ps) {
    pos(lam_nm, 'Wavelength'); pos(L_km, 'Fiber length');
    req(isFinite(D), 'Dispersion parameter D must be a number');
    req(isFinite(dlam_nm) && dlam_nm >= 0, 'Source spectral width must be zero or positive');
    req(isFinite(tfwhm_ps) && tfwhm_ps >= 0, 'Pulse width must be zero or positive');
    var lam = lam_nm * 1e-9, Dsi = D * 1e-6;                                      // s/m^2
    var b2 = -Dsi * lam * lam / (2 * PI * C);                                      // s^2/m
    var b2_ps2km = b2 * 1e24 * 1e3;
    var r = { beta2_ps2_km: b2_ps2km, total_ps_nm: D * L_km, regime: D > 0 ? 'anomalous (D > 0, β₂ < 0)' : D < 0 ? 'normal (D < 0, β₂ > 0)' : 'zero dispersion',
      spread_ps: null, LD_km: null, tout_ps: null, broadening: null, chirp_rate_ps_nm: null };
    if (dlam_nm > 0) r.spread_ps = Math.abs(D) * L_km * dlam_nm;
    if (tfwhm_ps > 0) {
      var T0 = tfwhm_ps / (2 * Math.sqrt(Math.LN2));
      if (b2_ps2km !== 0) { r.LD_km = T0 * T0 / Math.abs(b2_ps2km); var x = L_km / r.LD_km; r.broadening = Math.sqrt(1 + x * x); }
      else { r.LD_km = Infinity; r.broadening = 1; }
      r.tout_ps = tfwhm_ps * r.broadening;
    }
    return r;
  }

  // 27. uniform fiber Bragg grating (Erdogan, J. Lightwave Technol. 15, 1277 (1997)). Bragg wavelength nm, n_eff,
  // index modulation amplitude dn (ac), mode overlap eta (0..1], grating length mm. kappa = pi*eta*dn/lambda_B;
  // R_peak = tanh^2(kappa L); zero-to-zero bandwidth = lambda_B^2/(pi n_eff L) * sqrt((kappa L)^2 + pi^2).
  // Sensitivities for silica fiber: dlam/dT = lambda_B*(alpha + xi), alpha = 0.55e-6/K, xi = 8.6e-6/K (Othonos, Rev. Sci. Instrum. 68, 4309 (1997), Ge-doped silica);
  // dlam/eps = lambda_B*(1 - p_e), p_e = 0.22.
  function fbg(lam_nm, neff, dn, eta, L_mm) {
    pos(lam_nm, 'Bragg wavelength'); pos(neff, 'Effective index'); pos(dn, 'Index modulation'); pos(L_mm, 'Grating length');
    req(isFinite(eta) && eta > 0 && eta <= 1, 'Mode overlap must be between 0 and 1');
    req(neff >= 1, 'Effective index must be at least 1');
    var lam = lam_nm * 1e-9, L = L_mm * 1e-3;
    var kappa = PI * eta * dn / lam, kL = kappa * L, t = Math.tanh(kL), R = t * t;
    var bw0 = lam * lam / (PI * neff * L) * Math.sqrt(kL * kL + PI * PI);
    return { period_nm: lam_nm / (2 * neff), kappa_per_mm: kappa * 1e-3, kL: kL, R: R, R_pct: R * 100,
      R_dB: R < 1 ? 10 * Math.log10(1 - R) : -Infinity, bw0_nm: bw0 * 1e9,
      dT_pm_K: lam_nm * 1e3 * (0.55e-6 + 8.6e-6), dEps_pm_ue: lam_nm * 1e3 * (1 - 0.22) * 1e-6 };
  }

  // 28. symmetric directional coupler (coupled-mode / supermode picture). lambda nm, supermode indices n_even > n_odd,
  // coupler length um, target cross-port fraction (0 = skip). Beat (full-transfer) length L_c = lambda/(2*dn);
  // kappa = pi*dn/lambda = pi/(2 L_c); P_cross = sin^2(kappa L), P_bar = cos^2(kappa L) (lossless, phase-matched).
  function dirCoupler(lam_nm, n_even, n_odd, L_um, target) {
    pos(lam_nm, 'Wavelength'); pos(n_even, 'Even supermode index'); pos(n_odd, 'Odd supermode index'); pos(L_um, 'Coupler length');
    req(n_even > n_odd, 'The even (symmetric) supermode index must be larger than the odd one');
    req(isFinite(target) && target >= 0 && target <= 1, 'Target cross fraction must be between 0 and 1');
    var lam_um = lam_nm * 1e-3, dn = n_even - n_odd, Lc = lam_um / (2 * dn), kappa = PI * dn / lam_um, kL = kappa * L_um;
    var cross = Math.pow(Math.sin(kL), 2);
    return { dn: dn, Lc_um: Lc, kappa_per_um: kappa, kL: kL, cross: cross, cross_pct: cross * 100, bar_pct: (1 - cross) * 100,
      cross_dB: cross > 0 ? 10 * Math.log10(cross) : -Infinity, bar_dB: cross < 1 ? 10 * Math.log10(1 - cross) : -Infinity,
      L50_um: Lc / 2, Ltarget_um: target > 0 ? Math.asin(Math.sqrt(target)) / kappa : null };
  }

  // 29. unbalanced Mach-Zehnder interferometer (Chrostowski & Hochberg 2015, ch. 4.5). lambda nm, group index n_g,
  // arm length difference dL um, effective (phase) index n_eff, target FSR (0 = skip) in nm or GHz.
  // FSR = lambda^2/(n_g dL) = c/(n_g dL) in frequency; dphi = 2 pi n_eff dL / lambda; ideal 50:50 splitter/combiner, lossless arms:
  // port A = cos^2(dphi/2), port B = sin^2(dphi/2). Nearest port-A peak from linearised phase (slope set by n_g).
  function mzi(lam_nm, ng, dL_um, neff, target, tunit) {
    pos(lam_nm, 'Wavelength'); pos(ng, 'Group index'); pos(dL_um, 'Path length difference'); pos(neff, 'Effective index');
    req(isFinite(target) && target >= 0, 'Target FSR must be zero (skip) or positive');
    req(tunit === 'nm' || tunit === 'GHz', 'Unknown FSR unit');
    var lam = lam_nm * 1e-9, dL = dL_um * 1e-6;
    var fsr_m = lam * lam / (ng * dL), fsr_Hz = C / (ng * dL), tau = ng * dL / C;
    var phi = 2 * PI * neff * dL / lam, orders = neff * dL / lam, wrapped = phi - 2 * PI * Math.floor(phi / (2 * PI));
    var A = Math.pow(Math.cos(phi / 2), 2), d = wrapped > PI ? wrapped - 2 * PI : wrapped;
    var dL_t = null;
    if (target > 0) dL_t = tunit === 'nm' ? lam * lam / (ng * target * 1e-9) * 1e6 : C / (ng * target * 1e9) * 1e6;
    return { fsr_nm: fsr_m * 1e9, fsr_GHz: fsr_Hz / 1e9, delay_ps: tau * 1e12, phi_rad: wrapped, phi_pi: wrapped / PI, orders: orders,
      portA: A, portB: 1 - A, portA_pct: A * 100, portB_pct: (1 - A) * 100,
      peakA_nm: lam_nm + d / (2 * PI) * fsr_m * 1e9, dL_target_um: dL_t };
  }

  // 30. fiber optic link loss budget (FOA / TIA-568 style). Powers in dBm, losses in dB. Total passive loss =
  // L*alpha + nc*lc + ns*ls (+ other); available budget = P_tx - S_rx; margin = budget - loss - safety margin.
  // Maximum reach = (budget - connectors - splices - other - safety margin) / alpha (0 if the fixed losses already exceed it).
  function linkBudget(ptx, srx, L_km, alpha, nc, lc, ns, ls, other, safety) {
    req(isFinite(ptx) && isFinite(srx), 'Transmitter power and receiver sensitivity must be numbers');
    req(ptx > srx, 'Transmitter power must be higher than the receiver sensitivity');
    req(isFinite(L_km) && L_km >= 0, 'Fiber length must be zero or positive'); pos(alpha, 'Fiber attenuation');
    [nc, ns].forEach(function (n) { req(isFinite(n) && n >= 0 && Math.floor(n) === n, 'Connector and splice counts must be whole numbers'); });
    [lc, ls, other, safety].forEach(function (x) { req(isFinite(x) && x >= 0, 'Losses and margins must be zero or positive'); });
    var fiber = L_km * alpha, conn = nc * lc, spl = ns * ls, loss = fiber + conn + spl + other, budget = ptx - srx;
    var margin = budget - loss - safety, fixed = conn + spl + other + safety;
    return { fiber_dB: fiber, conn_dB: conn, splice_dB: spl, loss_dB: loss, budget_dB: budget, prx_dBm: ptx - loss,
      margin_dB: margin, ok: margin >= 0, maxL_km: Math.max(0, (budget - fixed) / alpha) };
  }

  // 31. OSNR of a chain of identical amplified spans (Agrawal 2010 ch. 7; Ramaswami et al. 2010). Each amplifier's gain G equals
  // the span loss; its ASE (both polarizations) in the reference bandwidth B_ref is P_ASE = NF*h*nu*G*B_ref (high-gain approximation).
  // N amplifiers add their ASE: OSNR = P_ch / (N*NF*h*nu*G*B_ref), i.e. OSNR[dB] = P_ch - L_span - NF - 10log10(N) - 10log10(h nu B_ref / 1 mW).
  // At 1550 nm and B_ref = 0.1 nm the last term is +57.96 dB (the "58 dB rule"). B_ref[Hz] = c*dlam/lam^2.
  function osnr(Pch_dBm, Lspan_dB, NF_dB, N, lam_nm, bref_nm, req_dB, nch) {
    req(isFinite(Pch_dBm), 'Launch power per channel must be a number'); pos(Lspan_dB, 'Span loss');
    req(isFinite(NF_dB) && NF_dB >= 0, 'Noise figure must be zero or positive (dB)');
    [N, nch].forEach(function (n) { req(isFinite(n) && n >= 1 && Math.floor(n) === n, 'Number of spans and channels must be whole numbers of at least 1'); });
    pos(lam_nm, 'Wavelength'); pos(bref_nm, 'Reference bandwidth');
    req(isFinite(req_dB) && req_dB >= 0, 'Required OSNR must be zero (skip) or positive');
    var lam = lam_nm * 1e-9, nu = C / lam, Bref = C * bref_nm * 1e-9 / (lam * lam);
    var k_dB = -10 * Math.log10(H * nu * Bref / 1e-3);
    var one = k_dB + Pch_dBm - Lspan_dB - NF_dB, tot = one - 10 * Math.log10(N);
    var ase1 = NF_dB + Lspan_dB - k_dB, nmax = null, margin = null;
    if (req_dB > 0) { margin = tot - req_dB; nmax = one >= req_dB ? Math.floor(Math.pow(10, (one - req_dB) / 10) + 1e-9) : 0; }
    return { const_dB: k_dB, bref_GHz: Bref / 1e9, nu_THz: nu / 1e12, osnr1_dB: one, osnr_dB: tot, ase1_dBm: ase1,
      aseN_dBm: ase1 + 10 * Math.log10(N), total_dBm: Pch_dBm + 10 * Math.log10(nch), margin_dB: margin,
      ok: margin === null ? null : margin >= 0, nmax: nmax };
  }

  // 32. Airy disk / diffraction-limited resolution of a uniformly illuminated circular aperture (Airy 1835; Born & Wolf,
  // Principles of Optics, 7th ed., sec. 8.5.2 and 8.6.2; Hecht, Optics, 5th ed., sec. 10.2.5). Intensity I(v) = I0 [2 J1(v)/v]^2 with
  // v = pi*D*r/(lambda*f) = 2*pi*NA*r/lambda. Dark rings at the zeros of J1: v = 3.8317, 7.0156 -> r = (v/pi)*lambda*N
  // (N = f/D = 1/(2NA)), i.e. r_Airy = 1.21967*lambda*N = 0.60983*lambda/NA (Rayleigh criterion). Half maximum at v = 1.61634 ->
  // FWHM = 1.02899*lambda*N. Abbe limit d = lambda/(2NA). Angular radius theta = 1.21967*lambda/D. Incoherent MTF cutoff
  // nu_c = 1/(lambda*N) = 2NA/lambda (Goodman, Introduction to Fourier Optics, 3rd ed., sec. 6.3), Nyquist pixel pitch = lambda*N/2.
  // Paraxial: NA = D/(2f) when NA is not given. 83.8 % of the energy lies inside the first dark ring.
  var J1Z1 = 3.8317059702075125, J1Z2 = 7.015586669815619, AIRY_HM = 1.6163399561827614;
  function airy(lam_nm, D_mm, f_mm, NA_in, pix_um) {
    pos(lam_nm, 'Wavelength'); pos(D_mm, 'Aperture diameter');
    req(isFinite(f_mm) && f_mm >= 0, 'Focal length must be zero (use NA) or positive');
    req(isFinite(NA_in) && NA_in >= 0 && NA_in <= 1.7, 'NA must be zero (use D and f) or between 0 and 1.7');
    req(f_mm > 0 || NA_in > 0, 'Give either a focal length or a numerical aperture');
    req(isFinite(pix_um) && pix_um >= 0, 'Pixel pitch must be zero (skip) or positive');
    var lam_um = lam_nm / 1000, NA = NA_in > 0 ? NA_in : D_mm / (2 * f_mm), N = 1 / (2 * NA);
    var r = (J1Z1 / PI) * lam_um * N, theta = (J1Z1 / PI) * lam_nm * 1e-9 / (D_mm * 1e-3);
    var nyq = lam_um * N / 2, out = { NA: NA, fnum: N, r_um: r, d_um: 2 * r, r2_um: (J1Z2 / PI) * lam_um * N,
      fwhm_um: (2 * AIRY_HM / PI) * lam_um * N, abbe_um: lam_um / (2 * NA), theta_urad: theta * 1e6, theta_arcsec: theta * R2D * 3600,
      cutoff_lpmm: 1000 / (lam_um * N), nyq_um: nyq, px_per_d: null, sampled: null, source: NA_in > 0 ? 'NA' : 'D/f' };
    if (pix_um > 0) { out.px_per_d = 2 * r / pix_um; out.sampled = pix_um <= nyq * (1 + 1e-12); }
    return out;
  }

  // 33. Single-layer anti-reflection coating at normal incidence (Hecht, Optics, 5th ed., sec. 9.7.1; Born & Wolf, Principles of
  // Optics, 7th ed., sec. 1.6.4; Macleod, Thin-Film Optical Filters, 4th ed., ch. 3). Interface amplitude coefficients
  // r1 = (n0 - n1)/(n0 + n1), r2 = (n1 - ns)/(n1 + ns); phase thickness delta = 2*pi*n1*d/lambda; Airy sum
  // R = (r1^2 + r2^2 + 2 r1 r2 cos 2delta) / (1 + r1^2 r2^2 + 2 r1 r2 cos 2delta). Quarter-wave layer d = m*lambda0/(4 n1), m odd,
  // gives R = ((n0 ns - n1^2)/(n0 ns + n1^2))^2, zero when n1 = sqrt(n0 ns). Window with two identical surfaces, no absorption,
  // incoherent multiple reflections: T = (1 - R)/(1 + R). Ignores dispersion of the indices and absorption.
  function arCoat(lam0_nm, n1, ns, n0, lam_nm, m) {
    pos(lam0_nm, 'Design wavelength'); pos(n1, 'Coating index'); pos(ns, 'Substrate index'); pos(n0, 'Ambient index');
    req(isFinite(lam_nm) && lam_nm >= 0, 'Evaluation wavelength must be zero (use design) or positive');
    req(isFinite(m) && m >= 1 && Math.round(m) === m && m % 2 === 1, 'Quarter-wave order must be an odd whole number (1, 3, 5, ...)');
    var lam = lam_nm > 0 ? lam_nm : lam0_nm, d = m * lam0_nm / (4 * n1);
    var r1 = (n0 - n1) / (n0 + n1), r2 = (n1 - ns) / (n1 + ns), c2 = Math.cos(4 * PI * n1 * d / lam);
    var R = (r1 * r1 + r2 * r2 + 2 * r1 * r2 * c2) / (1 + r1 * r1 * r2 * r2 + 2 * r1 * r2 * c2);
    var q = (n0 * ns - n1 * n1) / (n0 * ns + n1 * n1), Rq = q * q, Rb = Math.pow((n0 - ns) / (n0 + ns), 2);
    return { ideal_n: Math.sqrt(n0 * ns), d_nm: d, opt_nm: n1 * d, lam_nm: lam, R_bare_pct: Rb * 100, R_design_pct: Rq * 100, R_pct: R * 100,
      gain_pct: (Rb - R) * 100, T_win_bare_pct: (1 - Rb) / (1 + Rb) * 100, T_win_pct: (1 - R) / (1 + R) * 100,
      reduction: R > 0 ? Rb / R : Infinity };
  }

  var PC = { photon: photon, power: power, dbRatio: dbRatio, gaussian: gaussian, thinLens: thinLens, numAperture: numAperture, vNumber: vNumber, snell: snell, pulse: pulse, neMgLN: neMgLN, ppln: ppln, fabryPerot: fabryPerot, grating: grating, lc: lc, qpm3: qpm3, pplnTemp: pplnTemp, ring: ring, sellmeier: sellmeier, coupling: coupling, gratingCoupler: gratingCoupler, tbp: tbp, photodiode: photodiode, focusSpot: focusSpot, cohLen: cohLen, neCLN: neCLN, slabN: slabN, tfln: tfln, eoVpi: eoVpi, coherence: coherence, dispersion: dispersion, fbg: fbg, dirCoupler: dirCoupler, mzi: mzi, linkBudget: linkBudget, osnr: osnr, airy: airy, arCoat: arCoat, consts: { C: C, H: H, QE: QE } };
  if (typeof module !== 'undefined' && module.exports) module.exports = PC; else root.PC = PC;
})(this);
