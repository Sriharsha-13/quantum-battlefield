// Quantum Plants vs. Zombies: Adventure 3, "Ground-state hunters" (2 levels, 3 qubits)
// a. H = X0X1 + X1X2 (the brief's X1X2 + X2X3 with Qiskit's 0-based labels). Tune RY angles by hand; one
//    Trotter step e^{-iφH/2} (two RXX gates) is planted too, and it can never change the energy: [H, H] = 0.
// b. H = X0X1 + X1X2 + ½(Z0 + Z1 + Z2). Adiabatic state preparation from |111⟩ along
//    H(s) = s·(X0X1 + X1X2) + (1 − s/2)·ΣZ, Trotterized with n steps of 1st or 2nd (Suzuki) order.
// The boss falls only to a shot-based energy estimate below its shield, measured as in VQE: X-basis shots
// for the X·X terms, Z-basis shots for the Z terms. Bit strings read q2 q1 q0, as in Qiskit.
(() => {
  "use strict";

  // ==QUANTUM CORE==
  const NQ = 3, DIM = 8;
  const SQ = Math.SQRT1_2;
  const C = (re, im = 0) => [re, im];
  const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const cabs2 = a => a[0] * a[0] + a[1] * a[1];
  const bitOf = (i, q) => (i >> q) & 1;
  const basisState = i => { const a = Array.from({ length: DIM }, () => C(0)); a[i] = C(1); return a; };
  // P is a Pauli string as {qubit: "X" | "Y" | "Z"}.
  function pauliApply(a, P) {
    const b = Array.from({ length: DIM }, () => C(0));
    for (let i = 0; i < DIM; i++) {
      let j = i, ph = C(1);
      for (const k in P) {
        const q = +k, bq = bitOf(i, q);
        if (P[k] === "X") j ^= 1 << q;
        else if (P[k] === "Z") { if (bq) ph = [-ph[0], -ph[1]]; }
        else { j ^= 1 << q; ph = cmul(ph, bq ? C(0, -1) : C(0, 1)); }
      }
      b[j] = cadd(b[j], cmul(ph, a[i]));
    }
    return b;
  }
  // exp(−i·angle/2·P): RX/RY/RZ for one-qubit P, RXX for P = X⊗X.
  function pauliRot(a, P, angle) {
    const pa = pauliApply(a, P), c = Math.cos(angle / 2), s = Math.sin(angle / 2);
    return a.map((v, i) => cadd([v[0] * c, v[1] * c], cmul(C(0, -s), pa[i])));
  }
  const RY = (a, q, t) => pauliRot(a, { [q]: "Y" }, t);
  const RZ = (a, q, t) => pauliRot(a, { [q]: "Z" }, t);
  const RXX = (a, q1, q2, t) => pauliRot(a, { [q1]: "X", [q2]: "X" }, t);
  const Xg = (a, q) => pauliApply(a, { [q]: "X" });
  const Hg = (a, q) => {
    const b = a.map(v => v.slice());
    for (let i = 0; i < DIM; i++) {
      if (bitOf(i, q)) continue;
      const j = i | (1 << q);
      b[i] = [(a[i][0] + a[j][0]) * SQ, (a[i][1] + a[j][1]) * SQ];
      b[j] = [(a[i][0] - a[j][0]) * SQ, (a[i][1] - a[j][1]) * SQ];
    }
    return b;
  };
  const expect = (a, P) => { const pa = pauliApply(a, P); return a.reduce((s, v, i) => s + v[0] * pa[i][0] + v[1] * pa[i][1], 0); };
  const energy = (a, H) => H.reduce((s, t) => s + t.c * expect(a, t.P), 0);
  const probs = a => a.map(cabs2);
  function sampleIndex(p) {
    let r = Math.random() * p.reduce((s, v) => s + v, 0);
    for (let i = 0; i < p.length; i++) { r -= p[i]; if (r < 0) return i; }
    return p.length - 1;
  }
  // Every term is all-X or all-Z, so two measurement settings cover the whole Hamiltonian.
  const groupOf = t => Object.values(t.P)[0];
  const rotateTo = (a, basis) => { let b = a; if (basis === "X") for (let q = 0; q < NQ; q++) b = Hg(b, q); return b; };
  // Energy carried by one shot k measured in `basis`: Σ c·Π(±1) over that group's terms.
  const shotEnergy = (k, basis, H) => H.filter(t => groupOf(t) === basis)
    .reduce((s, t) => s + t.c * Object.keys(t.P).reduce((pr, q) => pr * (1 - 2 * bitOf(k, +q)), 1), 0);
  // Real symmetric eigensolver (Jacobi); H here has only X and Z terms, so its matrix is real.
  function spectrum(H) {
    const A = [];
    for (let k = 0; k < DIM; k++) {
      let col = Array.from({ length: DIM }, () => C(0));
      for (const t of H) { const p = pauliApply(basisState(k), t.P); col = col.map((v, i) => cadd(v, [t.c * p[i][0], t.c * p[i][1]])); }
      A.push(col.map(v => v[0]));
    }
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let i = 0; i < DIM; i++) for (let j = i + 1; j < DIM; j++) off += A[i][j] ** 2;
      if (off < 1e-22) break;
      for (let p = 0; p < DIM; p++) for (let q = p + 1; q < DIM; q++) {
        if (Math.abs(A[p][q]) < 1e-15) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < DIM; k++) { const x = A[k][p], y = A[k][q]; A[k][p] = c * x - s * y; A[k][q] = s * x + c * y; }
        for (let k = 0; k < DIM; k++) { const x = A[p][k], y = A[q][k]; A[p][k] = c * x - s * y; A[q][k] = s * x + c * y; }
      }
    }
    return A.map((r, i) => r[i]).sort((x, y) => x - y);
  }
  // Adiabatic path for level b: H(s) = s·(X0X1 + X1X2) + (1 − s/2)·ΣZ, from |111⟩ (ground state at s = 0).
  function trotterStep(a, dt, s, order) {
    const cz = 1 - 0.5 * s, cx = s;
    const zLayer = (b, f) => { for (let q = 0; q < NQ; q++) b = RZ(b, q, 2 * cz * dt * f); return b; };
    let b = a;
    b = zLayer(b, order === 2 ? 0.5 : 1);
    b = RXX(b, 0, 1, 2 * cx * dt);
    b = RXX(b, 1, 2, 2 * cx * dt);
    if (order === 2) b = zLayer(b, 0.5);
    return b;
  }
  function adiabatic(a, T, n, order) {
    let b = a;
    for (let k = 0; k < n; k++) b = trotterStep(b, T / n, (k + 0.5) / n, order);
    return b;
  }
  // ==END CORE==

  const PI = Math.PI;
  const XX = [{ c: 1, P: { 0: "X", 1: "X" } }, { c: 1, P: { 1: "X", 2: "X" } }];
  const HFIELD = [...XX, { c: 0.5, P: { 0: "Z" } }, { c: 0.5, P: { 1: "Z" } }, { c: 0.5, P: { 2: "Z" } }];
  const LEVELS = [
    { key: "a", name: "Commuting ladder", H: XX, hText: "H = X₀X₁ + X₁X₂", bar: -1.8, groups: ["X"], volley: 6,
      params: [["θ0", "RY on q0"], ["θ1", "RY on q1"], ["θ2", "RY on q2"], ["φ", "Trotter step RXX(φ) = e^(−iφH/2)"]],
      intro: "Find the lowest-energy state of H = X₀X₁ + X₁X₂ by tuning the RY angles. Each X·X term is −1 when its two qubits point opposite ways along x. The planted Trotter step e^(−iφH/2) never changes the energy: H commutes with itself." },
    { key: "b", name: "Trotter–Suzuki", H: HFIELD, hText: "H = X₀X₁ + X₁X₂ + ½(Z₀ + Z₁ + Z₂)", bar: -2.2, groups: ["X", "Z"], volley: 200,
      intro: "A transverse field makes the terms non-commuting. Start in |111⟩, the ground state of the field, and slowly switch on the X·X couplings (adiabatic theorem). Choose the total time T, the number of Trotter steps n and the Suzuki order." }
  ];
  const SHOTS = 8000;
  const volley = () => lvl().volley;   // shots per register: few for level a, a larger job for the noisier b and c

  const $ = id => document.getElementById(id);
  const canvas = $("lawn"), ctx = canvas.getContext("2d");
  const specCv = $("specCv"), sctx = specCv.getContext("2d");
  const landCv = $("landCv"), lctx = landCv.getContext("2d");

  const W = 1200, H = 470, TOP = 44, LANE_H = 116, DATA_H = 3 * LANE_H, BOT = TOP + DATA_H;
  const X0 = 100, XM = 584, XE = 1172;
  const FONT = '"Lilita One", "Arial Rounded MT Bold", "Trebuchet MS", system-ui, sans-serif';
  const LANE_END = 15, MEASURE_X = 5, QSPEED = 1, EMIT = 2, BOSS_SPEED = 0.2, BOSS_LAB_SPEED = 0.12, RAD = 0.3;
  const xPos = u => (u <= MEASURE_X ? X0 + (u * (XM - X0)) / MEASURE_X : XM + ((u - MEASURE_X) * (XE - XM)) / (LANE_END - MEASURE_X));
  const laneY = q => TOP + q * LANE_H + LANE_H / 2;
  const lanes = [0, 1, 2];

  const G = {
    level: 0, p: {}, T: 2, n: 2, order: 1, sel: null, preview: false,
    phase: "plan", baseHP: 10, gateHP: 10, score: 0, regs: [], boss: null, splats: [], floats: [],
    emitT: 0, groupIdx: 0, paused: false, nextId: 1, gateFlash: 0, houseFlash: 0, monT: 0,
    stream: {}, duel: null, spec: [], cols: [], lastShot: ""
  };
  const lvl = () => LEVELS[G.level];
  const gateAlive = () => G.gateHP > 0;

  // ---------- Formatting ----------
  const gcd = (a, b) => (b ? gcd(b, a % b) : Math.abs(a));
  function fmtAng(t) {
    const k = Math.round(t / (PI / 12));
    if (Math.abs(t - (k * PI) / 12) > 0.004) return t.toFixed(2);
    if (k === 0) return "0";
    const g = gcd(k, 12), p = k / g, q = 12 / g;
    const num = `${p < 0 ? "−" : ""}${Math.abs(p) === 1 ? "" : Math.abs(p)}π`;
    return q === 1 ? num : `${num}/${q}`;
  }
  const fmtE = e => (e < 0 ? "−" : "") + Math.abs(e).toFixed(3);
  const sub = s => s.replace(/[0-9]/g, d => "₀₁₂₃₄₅₆₇₈₉"[d]);

  // ---------- Circuits ----------
  // Each column: { u, tiles: [{ qs, label, val, param, locked }], op(a) }.
  function buildCols() {
    const L = lvl(), P = G.p;
    if (L.key === "a") {
      return [
        { u: 1.0, tiles: lanes.map(q => ({ qs: [q], label: "RY", val: fmtAng(P[`θ${q}`]), param: `θ${q}` })),
          op: a => lanes.reduce((b, q) => RY(b, q, P[`θ${q}`]), a) },
        { u: 2.3, tiles: [{ qs: [0, 1], label: "RXX", val: fmtAng(P.φ), param: "φ" }], op: a => RXX(a, 0, 1, P.φ) },
        { u: 3.3, tiles: [{ qs: [1, 2], label: "RXX", val: fmtAng(P.φ), param: "φ" }], op: a => RXX(a, 1, 2, P.φ) }
      ];
    }
    // Level b: one Trotter step is drawn; the register receives all n steps as it leaves the block.
    const init = { u: 0.65, tiles: lanes.map(q => ({ qs: [q], label: "X", locked: true })), op: a => lanes.reduce((b, q) => Xg(b, q), a) };
    const z = f => ({ tiles: lanes.map(q => ({ qs: [q], label: "RZ", val: f, locked: true })), op: a => a });
    const cols = [init, { u: 1.5, ...z(G.order === 2 ? "½·step" : "step") },
      { u: 2.3, tiles: [{ qs: [0, 1], label: "RXX", val: "step", locked: true }], op: a => a },
      { u: 3.1, tiles: [{ qs: [1, 2], label: "RXX", val: "step", locked: true }], op: a => a }];
    if (G.order === 2) cols.push({ u: 3.9, ...z("½·step") });
    cols[cols.length - 1].op = a => adiabatic(a, G.T, G.n, G.order);
    cols[cols.length - 1].blockEnd = true;
    return cols;
  }
  const prepare = () => G.cols.reduce((a, c) => c.op(a), basisState(0));
  const exactE = () => energy(prepare(), lvl().H);
  // Exact (Trotter-error-free) adiabatic result for level b, for comparison.
  const exactAdiabatic = T => energy(adiabatic(basisState(7), T, 240, 2), lvl().H);

  function qiskitCode() {
    const L = lvl();
    const lab = t => { const s = ["I", "I", "I"]; for (const q in t.P) s[2 - q] = t.P[q]; return s.join(""); };
    const lines = ["import { QuantumCircuit, SparsePauliOp } from 'qiskit';", "",
      `// ${L.hText}; Pauli labels read q2 q1 q0`,
      `const H = SparsePauliOp.from_list([${L.H.map(t => `["${lab(t)}", ${t.c}]`).join(", ")}]);`,
      "const qc = QuantumCircuit(3);"];
    const P = G.p, f = x => x.toFixed(4);
    if (L.key === "a") {
      lanes.forEach(q => lines.push(`qc.ry(${f(P[`θ${q}`])}, ${q});`));
      lines.push(`qc.rxx(${f(P.φ)}, 0, 1); qc.rxx(${f(P.φ)}, 1, 2);   // one Trotter step of e^(-i φ H / 2)`);
    } else {
      lines.push("qc.x(0); qc.x(1); qc.x(2);   // |111⟩, ground state of the field term",
        `const T = ${G.T}, n = ${G.n}, dt = T / n;   // ${G.order === 2 ? "2nd-order (Suzuki)" : "1st-order (Lie–Trotter)"} steps`,
        "for (let k = 0; k < n; k++) {",
        "  const s = (k + 0.5) / n, cz = 1 - 0.5 * s, cx = s;",
        `  for (const q of [0, 1, 2]) qc.rz(2 * cz * dt${G.order === 2 ? " / 2" : ""}, q);`,
        "  qc.rxx(2 * cx * dt, 0, 1); qc.rxx(2 * cx * dt, 1, 2);",
        ...(G.order === 2 ? ["  for (const q of [0, 1, 2]) qc.rz(2 * cz * dt / 2, q);"] : []),
        "}");
    }
    lines.push("", "// Energy: estimate <H> on the Qollab backend (Estimator), or measure two groups:",
      "// X-basis shots (H on every qubit, then measure) for the X·X terms, Z-basis shots for the Z terms.");
    return lines.join("\n");
  }

  // ---------- Levels ----------
  function fitCanvas() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = W * dpr; canvas.height = H * dpr; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    specCv.width = 420 * dpr; specCv.height = 230 * dpr; sctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    landCv.width = 420 * dpr; landCv.height = 170 * dpr; lctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }
  function say(t) { $("msg").textContent = t; }

  function loadLevel(i, keepStats) {
    G.level = i;
    const L = lvl();
    G.p = {};
    (L.params || []).forEach(([nm]) => { G.p[nm] = 0; });
    G.sel = L.params ? L.params[0][0] : "T";
    G.T = 2; G.n = 2; G.order = 1;
    if (!keepStats) { G.baseHP = 10; G.score = 0; }
    G.gateHP = 10; G.paused = false; G.splats = []; G.floats = []; G.regs = []; G.groupIdx = 0; G.emitT = 0;
    G.duel = null;
    G.spec = spectrum(L.H);
    G.boss = { x: 13.6, prevX: 13.6, flash: 0, hitGate: false };
    G.phase = "plan";
    buildParams();
    $("duelInfo").textContent = "No energy estimate yet.";
    say(`Level ${L.key}: ${L.name}. ${L.intro} Press Ready when your circuit is set.`);
    circuitChanged();
    autoTutorial();
  }
  function startWave() {
    if (G.phase !== "plan") return;
    G.phase = "battle"; G.emitT = EMIT;
    say(`Level ${lvl().key} is moving. Each register is a volley of ${volley()} shots, measured in the X basis${lvl().groups.length > 1 ? " or the Z basis, alternately" : ""}.`);
    updateHUD();
  }
  function resetStream() { G.stream = {}; lvl().groups.forEach(g => { G.stream[g] = { n: 0, s: 0, s2: 0 }; }); }
  // Estimate ⟨H⟩ = Σ over groups of the mean shot energy, with its standard error. The variance gets a
  // floor of (range²/4)/n, the largest variance a shot in that group can have, divided by n, so a short
  // lucky streak of identical shots cannot claim zero error. It is negligible for an 8000-shot duel.
  const groupRange = g => 2 * lvl().H.filter(t => groupOf(t) === g).reduce((s, t) => s + Math.abs(t.c), 0);
  function estimateFrom(stats) {
    let e = 0, v = 0, ok = true;
    for (const g of lvl().groups) {
      const st = stats[g];
      if (!st || st.n < 2) { ok = false; continue; }
      const m = st.s / st.n, R = groupRange(g);
      e += m; v += (Math.max(0, st.s2 / st.n - m * m) + (R * R) / 4 / st.n) / st.n;
    }
    return ok ? { e, se: Math.sqrt(v) } : null;
  }
  // The live stream defeats the boss once ⟨H⟩ + 3σ < shield (about 99.7% confidence) with enough shots per basis.
  // It re-tests after every shot, so a looser 2σ rule lets lucky streaks through for states just above the shield.
  const MIN_SHOTS = 10;
  function streamProof() {
    const est = estimateFrom(G.stream), n = Math.min(...lvl().groups.map(g => (G.stream[g] ? G.stream[g].n : 0)));
    if (!est) return { est: null, n, done: false };
    return { ...est, est, n, done: n >= MIN_SHOTS && est.e + 3 * est.se < lvl().bar };
  }

  // ---------- Entities ----------
  function emit() {
    const basis = lvl().groups[G.groupIdx++ % lvl().groups.length];
    G.regs.push({ id: G.nextId++, x: 0.25, prevX: 0.25, amps: basisState(0), next: 0, measured: false, basis, k: null, e: null, flash: 0 });
  }
  function passCols(r) {
    while (r.next < G.cols.length && r.x >= G.cols[r.next].u) { r.amps = G.cols[r.next].op(r.amps); r.next++; }
  }
  // Each register is a volley of volley() shots of the same circuit (like one small backend job).
  // The lanes show the first shot's outcome; r.e is the volley's average shot energy.
  function measureReg(r) {
    const p = probs(rotateTo(r.amps, r.basis)), st = G.stream[r.basis];
    let sum = 0;
    for (let i = 0; i < volley(); i++) {
      const k = sampleIndex(p), e = shotEnergy(k, r.basis, lvl().H);
      if (i === 0) r.k = k;
      sum += e; st.n++; st.s += e; st.s2 += e * e;
    }
    r.e = sum / volley();
    r.measured = true; r.flash = 0.3; G.gateFlash = 0.25;
    G.lastShot = `${r.basis} basis, ${volley()} shots, first ${lanes.slice().reverse().map(q => sym(r, q)).join("")}, average shot energy ${fmtE(r.e)}`;
  }
  const sym = (r, q) => (r.basis === "Z" ? String(bitOf(r.k, q)) : (bitOf(r.k, q) ? "−" : "+"));
  function splat(u, q, color) { G.splats.push({ x: xPos(u), y: laneY(q) - 8, r: 6, life: 0.4, color }); }
  function damageGate(n, t) {
    if (!gateAlive()) return;
    G.gateHP = Math.max(0, G.gateHP - n); G.gateFlash = 0.4;
    say(G.gateHP === 0 ? "Measurement gate destroyed: no more shots, so no energy estimates until it is rebuilt next level." : `${t}. Gate HP ${G.gateHP}.`);
  }
  function damageBase(n, t) {
    G.baseHP = Math.max(0, G.baseHP - n); G.houseFlash = 0.8;
    say(`${t}. Base HP ${G.baseHP}.`);
    if (G.baseHP === 0) { G.phase = "over"; say("Base lost. Restart to play again."); }
  }
  function defeatBoss(t) {
    G.score += 50;
    lanes.forEach(q => splat(G.boss.x, q, "#ffd23f"));
    G.boss = null; say(t);
    G.score += 25; G.phase = "cleared";
    const from = G.level, last = from === LEVELS.length - 1;
    setTimeout(() => {
      if (G.phase !== "cleared" || G.level !== from) return;
      loadLevel(last ? G.level : G.level + 1, true);
      if (last) say("Adventure 3 complete: you prepared ground states by hand and by Trotterized adiabatic evolution.");
    }, 2200);
  }

  function update(dt) {
    G.gateFlash = Math.max(0, G.gateFlash - dt); G.houseFlash = Math.max(0, G.houseFlash - dt);
    for (const s of G.splats) { s.life -= dt; s.r += dt * 50; }
    G.splats = G.splats.filter(s => s.life > 0);
    for (const f of G.floats) { f.life -= dt; f.y -= dt * 30; }
    G.floats = G.floats.filter(f => f.life > 0);
    if (G.phase !== "battle" || G.paused) return;

    G.emitT += dt;
    while (G.emitT >= EMIT) { G.emitT -= EMIT; emit(); }
    for (const r of G.regs) {
      r.prevX = r.x; r.x += QSPEED * dt; r.flash = Math.max(0, r.flash - dt);
      passCols(r);
      if (!r.measured && gateAlive() && r.x >= MEASURE_X) measureReg(r);
      if (r.x > LANE_END + 0.5) r.dead = true;
    }
    const proof = streamProof();
    if (G.boss && proof.done) {
      defeatBoss(`Your live shots proved ⟨H⟩ ≈ ${fmtE(proof.e)} ± ${proof.se.toFixed(3)} sits below the shield at ${fmtE(lvl().bar)} (3σ, about 99.7% confidence). Boss defeated.`);
    }
    const B = G.boss;
    if (B) {
      B.prevX = B.x; B.x -= (B.x < MEASURE_X ? BOSS_LAB_SPEED : BOSS_SPEED) * dt; B.flash = Math.max(0, B.flash - dt);
      if (!B.hitGate && B.x <= MEASURE_X) { B.hitGate = true; if (gateAlive()) damageGate(2, "The Hamiltonian boss shoved through the measurement gate"); }
      if (B.x <= 0.3) {
        damageBase(3, "The Hamiltonian boss reached the base and regroups at the far end");
        Object.assign(B, { x: 13.6, prevX: 13.6, hitGate: false });
      }
      // Measured shots that reach the boss splash their energy on it; the running estimate is what counts.
      for (const r of G.regs) {
        if (!r.measured || r.dead || B.x < MEASURE_X) continue;
        const s = r.prevX - B.prevX, e = r.x - B.x;
        if (!(Math.abs(e) <= RAD || s * e < 0)) continue;
        r.dead = true; B.flash = 0.15;
        const good = r.e < lvl().bar / lvl().groups.length;
        lanes.forEach(q => splat(B.x, q, good ? "#b6f27a" : "#c9a7ff"));
        G.floats.push({ x: xPos(B.x) + 115, y: laneY(1), life: 1.1, good,
          s: lvl().groups.length > 1 ? `${r.basis} terms, avg ${fmtE(r.e)}` : `${volley()} shots, avg ${fmtE(r.e)}` });
      }
    }
    G.regs = G.regs.filter(r => !r.dead);
  }

  // ---------- Energy duel ----------
  function duel() {
    if (G.phase !== "battle" || !G.boss) return;
    if (!gateAlive()) { say("The measurement gate is destroyed, so no shots can be taken."); return; }
    const a = prepare(), L = lvl(), per = SHOTS / L.groups.length, stats = {}, termSum = {};
    for (const g of L.groups) {
      const p = probs(rotateTo(a, g)), st = { n: 0, s: 0, s2: 0 };
      for (let k = 0; k < per; k++) {
        const idx = sampleIndex(p), e = shotEnergy(idx, g, L.H);
        st.n++; st.s += e; st.s2 += e * e;
        L.H.filter(t => groupOf(t) === g).forEach((t, j) => {
          const key = `${g}${j}`, v = Object.keys(t.P).reduce((pr, q) => pr * (1 - 2 * bitOf(idx, +q)), 1);
          termSum[key] = termSum[key] || { t, s: 0, n: 0 }; termSum[key].s += v; termSum[key].n++;
        });
      }
      stats[g] = st;
    }
    const est = estimateFrom(stats);
    G.duel = est;
    const terms = Object.values(termSum).map(o => `⟨${Object.entries(o.t.P).map(([q, P]) => sub(`${P}${q}`)).join("")}⟩ ≈ ${fmtE(o.s / o.n)}`).join(", ");
    const head = `${SHOTS} shots (${L.groups.map(g => `${per} in ${g}`).join(", ")}): ⟨H⟩ ≈ ${fmtE(est.e)} ± ${est.se.toFixed(3)}. ${terms}.`;
    if (est.e < L.bar) {
      $("duelInfo").textContent = `${head} Below the shield at ${fmtE(L.bar)}: boss defeated.`;
      defeatBoss("Energy duel won: your state sits below the boss's shield.");
    } else {
      G.boss.x = Math.max(0.4, G.boss.x - 1.5);
      $("duelInfo").textContent = `${head} Not below ${fmtE(L.bar)}. The boss surges forward.`;
      say("Energy duel lost. Lower the energy and try again.");
    }
    drawSpectrum();
    updateHUD();
  }

  function circuitChanged() {
    G.cols = buildCols();
    resetStream();
    G.duel = null;
    $("qiskitCode").textContent = qiskitCode();
    syncParams();
    drawSpectrum();
    drawLandscape();
    updateHUD();
  }

  // ---------- Panels ----------
  function buildParams() {
    const L = lvl();
    let h = "";
    if (L.params) {
      h = L.params.map(([nm, desc]) => `<div class="prow" data-p="${nm}"><span class="pname">${nm}</span>`
        + `<input type="range" min="0" max="${(2 * PI).toFixed(3)}" step="0.005" data-p="${nm}" aria-label="${nm}: ${desc}">`
        + `<span class="pval" id="pv-${nm}"></span><span class="pdesc">${desc}</span></div>`).join("");
    } else {
      h = `<div class="prow" data-p="T"><span class="pname">T</span><input type="range" id="Tsl" min="0.5" max="10" step="0.5" aria-label="Total evolution time T"><span class="pval" id="pv-T"></span><span class="pdesc">total evolution time</span></div>`
        + `<div class="row"><span class="inline-label">Trotter steps n</span><div class="seg" id="nSel">${[1, 2, 4, 8, 12, 16, 24].map(n => `<button data-n="${n}">${n}</button>`).join("")}</div></div>`
        + `<div class="row"><span class="inline-label">Order</span><div class="seg" id="oSel"><button data-o="1">1st: Lie–Trotter</button><button data-o="2">2nd: Suzuki</button></div></div>`;
    }
    $("params").innerHTML = h;
    $("zeroP").classList.toggle("hidden", !L.params);
    $("hHint").textContent = L.key === "a" ? "Tune by hand. ⟨X⟩ of RY(θ)|0⟩ is sin θ, so E = sin θ₀ sin θ₁ + sin θ₁ sin θ₂ before the RXX step."
      : "All gate angles follow from T, n and the schedule s = (k + ½)/n: RZ angle 2(1 − s/2)Δt, RXX angle 2sΔt, with Δt = T/n.";
    document.querySelectorAll("#params input[data-p]").forEach(inp => inp.addEventListener("input", () => {
      G.p[inp.dataset.p] = +inp.value; G.sel = inp.dataset.p; circuitChanged();
    }));
    document.querySelectorAll("#params .prow").forEach(row => row.addEventListener("pointerdown", () => { G.sel = row.dataset.p; syncParams(); drawLandscape(); }));
    if (!L.params) {
      $("Tsl").addEventListener("input", e => { G.T = +e.target.value; circuitChanged(); });
      document.querySelectorAll("#nSel button").forEach(b => b.addEventListener("click", () => { G.n = +b.dataset.n; circuitChanged(); }));
      document.querySelectorAll("#oSel button").forEach(b => b.addEventListener("click", () => { G.order = +b.dataset.o; circuitChanged(); }));
    }
    document.querySelectorAll("#levelSel button").forEach((b, i) => b.classList.toggle("on", i === G.level));
  }
  function syncParams() {
    const L = lvl();
    if (L.params) {
      L.params.forEach(([nm]) => {
        const inp = document.querySelector(`#params input[data-p="${nm}"]`);
        if (inp && Math.abs(+inp.value - G.p[nm]) > 1e-9) inp.value = G.p[nm];
        $(`pv-${nm}`).textContent = fmtAng(G.p[nm]);
      });
    } else {
      $("Tsl").value = G.T; $("pv-T").textContent = G.T.toFixed(1);
      document.querySelectorAll("#nSel button").forEach(b => b.classList.toggle("on", +b.dataset.n === G.n));
      document.querySelectorAll("#oSel button").forEach(b => b.classList.toggle("on", +b.dataset.o === G.order));
    }
    document.querySelectorAll("#params .prow").forEach(r => r.classList.toggle("sel", r.dataset.p === G.sel));
  }
  function updateHUD() {
    const L = lvl();
    $("levelChip").textContent = `${L.key} · ${G.level + 1} of ${LEVELS.length}`;
    $("baseHP").textContent = `${G.baseHP} / 10`;
    $("gateHP").textContent = gateAlive() ? `${G.gateHP} / 10` : "destroyed";
    $("score").textContent = G.score;
    $("shieldChip").textContent = `⟨H⟩ < ${fmtE(L.bar)}`;
    const live = estimateFrom(G.stream);
    $("liveChip").textContent = live ? `${fmtE(live.e)} ± ${live.se.toFixed(2)}` : "no shots";
    $("hText").textContent = L.hText;
    $("pause").textContent = G.paused ? "Resume" : "Pause";
    $("ready").classList.toggle("hidden", G.phase !== "plan");
    $("ready").textContent = `Ready: start level ${L.key}`;
    $("duel").disabled = G.phase !== "battle" || !G.boss || !gateAlive();
    $("duel").textContent = `Energy duel: estimate ⟨H⟩ with ${SHOTS} shots`;
    $("liveInfo").textContent = live ? `Live shots since your last change: ${Object.entries(G.stream).map(([g, s]) => `${s.n} ${g}`).join(", ")}. ⟨H⟩ ≈ ${fmtE(live.e)} ± ${live.se.toFixed(2)}. ${G.lastShot ? `Last: ${G.lastShot}.` : ""}`
      : `Live estimate needs at least 2 shots in ${L.groups.length > 1 ? "each basis" : "the X basis"}. ${G.lastShot ? `Last: ${G.lastShot}.` : ""}`;
    const pv = G.preview;
    $("exactInfo").textContent = !pv ? "Simulator preview hidden. On real hardware you only get shots; turn it on to see exact energies."
      : L.key === "b" ? `Exact ⟨H⟩ of this Trotter circuit: ${fmtE(exactE())}. Same T with no Trotter error: ${fmtE(exactAdiabatic(G.T))}. Ground state: ${fmtE(G.spec[0])}.`
      : `Exact ⟨H⟩ of your circuit: ${fmtE(exactE())}. Ground state: ${fmtE(G.spec[0])}${Math.abs(G.spec[0] - G.spec[1]) < 1e-6 ? " (two-fold degenerate)" : ""}.`;
  }

  // ---------- Energy panels ----------
  const isDark = () => window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
  function drawSpectrum() {
    const c = sctx, w = 420, h = 230, L = lvl();
    c.clearRect(0, 0, w, h);
    const ink = isDark() ? "#f6e7c1" : "#3b2412", faint = isDark() ? "rgba(246,231,193,.45)" : "rgba(59,36,18,.35)";
    const lo = G.spec[0] - 0.4, hi = G.spec[DIM - 1] + 0.4, Y = e => 14 + ((hi - e) / (hi - lo)) * (h - 28);
    c.font = `12px ${FONT}`; c.textBaseline = "middle";
    let lastY = -99;
    G.spec.forEach((e, i) => {
      const y = Y(e);
      c.strokeStyle = i === 0 ? "#2f9e1f" : faint; c.lineWidth = i === 0 ? 3 : 1.5;
      c.beginPath(); c.moveTo(60, y); c.lineTo(250, y); c.stroke();
      if (Math.abs(y - lastY) > 12) { c.fillStyle = i === 0 ? "#2f9e1f" : ink; c.textAlign = "right"; c.fillText(fmtE(e), 54, y); lastY = y; }
    });
    c.fillStyle = ink; c.textAlign = "left"; c.fillText("energy levels of H", 60, 8);
    c.save(); c.setLineDash([6, 4]); c.strokeStyle = "#c8402c"; c.lineWidth = 2;
    c.beginPath(); c.moveTo(40, Y(L.bar)); c.lineTo(410, Y(L.bar)); c.stroke(); c.restore();
    c.fillStyle = "#c8402c"; c.textAlign = "left"; c.fillText(`shield ${fmtE(L.bar)}`, 150, Y(L.bar) - 9);
    const mark = (x, e, col, lab, se) => {
      const y = Y(e);
      if (se) { c.strokeStyle = col; c.lineWidth = 2; c.beginPath(); c.moveTo(x, Y(e - 2 * se)); c.lineTo(x, Y(e + 2 * se)); c.stroke(); }
      c.fillStyle = col; c.beginPath(); c.arc(x, y, 6, 0, 2 * PI); c.fill();
      c.textAlign = "left"; c.fillText(lab, x + 10, y);
    };
    if (G.preview) mark(280, exactE(), "#2f9e1f", "exact", 0);
    const live = estimateFrom(G.stream);
    if (live) mark(320, live.e, "#3a8fd0", "live", live.se);
    if (G.duel) mark(360, G.duel.e, "#e08a1e", "duel", G.duel.se);
  }
  function drawLandscape() {
    const c = lctx, w = 420, h = 170, L = lvl();
    c.clearRect(0, 0, w, h);
    $("landWrap").classList.toggle("hidden", !G.preview);
    if (!G.preview) return;
    const ink = isDark() ? "#f6e7c1" : "#3b2412";
    let xs, ys, ys2 = null, cur, label;
    if (L.params) {
      const nm = G.sel, keep = G.p[nm];
      xs = Array.from({ length: 121 }, (_, i) => (i / 120) * 2 * PI);
      ys = xs.map(x => { G.p[nm] = x; G.cols = buildCols(); return exactE(); });
      G.p[nm] = keep; G.cols = buildCols();
      cur = keep; label = `E vs ${nm} (others fixed), 0 to 2π`;
    } else {
      const keep = G.T;
      xs = Array.from({ length: 39 }, (_, i) => 0.5 + i * 0.25);
      ys = xs.map(x => { G.T = x; G.cols = buildCols(); return exactE(); });
      ys2 = xs.map(x => exactAdiabatic(x));
      G.T = keep; G.cols = buildCols();
      cur = keep; label = `E vs T with n = ${G.n}, order ${G.order} (dashed: no Trotter error)`;
    }
    const lo = Math.min(G.spec[0], ...ys) - 0.2, hi = Math.max(...ys, ...(ys2 || []), L.bar) + 0.2;
    const X = x => 40 + ((x - xs[0]) / (xs[xs.length - 1] - xs[0])) * (w - 56), Y = e => 20 + ((hi - e) / (hi - lo)) * (h - 36);
    c.save(); c.setLineDash([6, 4]); c.strokeStyle = "#c8402c"; c.lineWidth = 1.5;
    c.beginPath(); c.moveTo(40, Y(L.bar)); c.lineTo(w - 16, Y(L.bar)); c.stroke(); c.restore();
    c.strokeStyle = "rgba(47,158,31,.6)"; c.lineWidth = 1; c.beginPath(); c.moveTo(40, Y(G.spec[0])); c.lineTo(w - 16, Y(G.spec[0])); c.stroke();
    const line = (arr, col, dash) => {
      c.save(); if (dash) c.setLineDash([4, 4]); c.strokeStyle = col; c.lineWidth = 2.5; c.beginPath();
      arr.forEach((e, i) => (i ? c.lineTo(X(xs[i]), Y(e)) : c.moveTo(X(xs[i]), Y(e)))); c.stroke(); c.restore();
    };
    if (ys2) line(ys2, "#7a4fc0", true);
    line(ys, "#3a8fd0", false);
    c.fillStyle = "#e08a1e"; c.beginPath(); c.arc(X(cur), Y(exactE()), 5, 0, 2 * PI); c.fill();
    c.fillStyle = ink; c.font = `12px ${FONT}`; c.textAlign = "left"; c.textBaseline = "middle";
    c.fillText(label, 40, 9);
    c.textAlign = "right"; c.fillText(fmtE(hi - 0.2), 36, Y(hi - 0.2)); c.fillText(fmtE(lo + 0.2), 36, Y(lo + 0.2));
  }

  // ---------- Drawing ----------
  function rr(x, y, w, h, r) { ctx.beginPath(); ctx.roundRect(x, y, w, h, r); }
  function circle(x, y, r) { ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); }
  function text(s, x, y, size, color, align = "center") {
    ctx.font = `${size}px ${FONT}`; ctx.fillStyle = color; ctx.textAlign = align; ctx.textBaseline = "middle"; ctx.fillText(s, x, y);
  }
  function drawLawn() {
    ctx.fillStyle = "#3d8a2a"; ctx.fillRect(0, 0, W, TOP);
    ctx.fillStyle = "#4c9c33";
    for (let x = 0; x < W; x += 34) { circle(x + 17, TOP - 2, 18); ctx.fill(); }
    for (const q of lanes) for (let u = 0; u < LANE_END; u++) {
      ctx.fillStyle = (q + u) % 2 ? "#6cba43" : "#7cc952";
      ctx.fillRect(xPos(u), TOP + q * LANE_H, xPos(u + 1) - xPos(u) + 1, LANE_H);
    }
    ctx.fillStyle = "rgba(255,255,240,.12)"; ctx.fillRect(xPos(0.4), TOP, xPos(MEASURE_X) - xPos(0.4), DATA_H);
    ctx.fillStyle = "#b98a4e"; ctx.fillRect(0, BOT, W, H - BOT);
    for (const u of [0, 5, 10, 15]) text(`x = ${u}`, xPos(u), H - 12, 12, "#5b3a1a");
    text("ansatz: prepare the trial state", (xPos(0.4) + xPos(MEASURE_X)) / 2, BOT + 17, 12, "#5b3a1a");
    text(`each register = ${volley()} shots, measured in ${lvl().groups.join(" / ")}`, (xPos(MEASURE_X) + xPos(LANE_END)) / 2, BOT + 17, 12, "#5b3a1a");
  }
  function drawHouse() {
    ctx.fillStyle = "#efdfb8"; ctx.fillRect(0, 0, 95, BOT);
    ctx.strokeStyle = "#cdb98b"; ctx.lineWidth = 2;
    for (let y = 10; y < BOT; y += 18) { ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(95, y); ctx.stroke(); }
    ctx.fillStyle = "#8a5a2b"; ctx.fillRect(88, 0, 9, BOT);
    const dy = TOP + DATA_H / 2;
    rr(24, dy - 38, 44, 76, 6); ctx.fillStyle = "#7a4a22"; ctx.fill();
    circle(60, dy, 4); ctx.fillStyle = "#ffd23f"; ctx.fill();
    text(`Base ${G.baseHP}`, 46, BOT + 18, 15, "#fff4cf");
    if (G.houseFlash > 0) { ctx.fillStyle = `rgba(220,50,40,${G.houseFlash * 0.5})`; ctx.fillRect(0, 0, 97, BOT); }
  }
  function drawSprout(x, y, t) {
    ctx.save(); ctx.translate(x, y + 6 + Math.sin(t * 2 + y) * 1.2); ctx.scale(0.5, 0.5);
    ctx.strokeStyle = "#2f7d1f"; ctx.lineWidth = 6;
    ctx.beginPath(); ctx.moveTo(0, 38); ctx.quadraticCurveTo(-6, 18, 0, 0); ctx.stroke();
    ctx.fillStyle = "#73d044"; ctx.strokeStyle = "#2f7d1f"; ctx.lineWidth = 3;
    rr(12, -26, 22, 20, 7); ctx.fill(); ctx.stroke();
    circle(0, -16, 20); ctx.fill(); ctx.stroke();
    circle(4, -24, 6); ctx.fillStyle = "#fff"; ctx.fill();
    circle(6, -24, 3); ctx.fillStyle = "#111"; ctx.fill();
    ctx.restore();
  }
  const TILE = { RY: ["#d2f1ff", "#1f6f9a", "#124a6a"], RZ: ["#f3e3ff", "#6b3fa0", "#46206e"], RXX: ["#ffe9a6", "#9a6a10", "#7a4a00"], X: ["#ffb347", "#9a5b10", "#5b3a1a"] };
  function drawCircuit(t) {
    const L = lvl();
    for (const q of lanes) {
      ctx.strokeStyle = "rgba(60,40,15,.45)"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(xPos(0.3), laneY(q)); ctx.lineTo(xPos(MEASURE_X), laneY(q)); ctx.stroke();
      drawSprout(xPos(0.15), laneY(q) - 4, t);
      text(`q${q}`, xPos(0.15) - 2, laneY(q) + 34, 12, "#2a4a12");
    }
    for (const col of G.cols) {
      const cx = xPos(col.u);
      for (const tl of col.tiles) {
        const [fill, stroke, ink] = TILE[tl.label];
        const y1 = laneY(Math.min(...tl.qs)) - 22, y2 = laneY(Math.max(...tl.qs)) + 22, w = tl.label === "X" ? 30 : 48;
        const yy1 = tl.label === "X" ? laneY(tl.qs[0]) - 15 : y1, hh = tl.label === "X" ? 30 : y2 - y1;
        rr(cx - w / 2, yy1, w, hh, 7); ctx.fillStyle = fill; ctx.fill();
        ctx.lineWidth = tl.param && tl.param === G.sel ? 4 : 2.5; ctx.strokeStyle = tl.param && tl.param === G.sel ? "#e08a1e" : stroke; ctx.stroke();
        const ym = (y1 + y2) / 2;
        if (tl.val) { text(tl.label, cx, ym - 9, 14, ink); text(tl.val, cx, ym + 9, 12, ink); } else text(tl.label, cx, ym, 15, ink);
        if (tl.locked) { rr(cx + w / 2 - 12, yy1 - 6, 9, 8, 2); ctx.fillStyle = "#8a5a2b"; ctx.fill(); }
      }
    }
    const bracket = (c1, c2, s) => {
      const a = xPos(G.cols[c1].u) - 26, b = xPos(G.cols[c2].u) + 26, y = 30;
      ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(a, y + 6); ctx.quadraticCurveTo(a, y, a + 8, y); ctx.lineTo(b - 8, y); ctx.quadraticCurveTo(b, y, b, y + 6); ctx.stroke();
      ctx.font = `13px ${FONT}`;
      const tw = ctx.measureText(s).width + 14;
      rr((a + b) / 2 - tw / 2, y - 20, tw, 18, 5); ctx.fillStyle = "rgba(255,248,220,.94)"; ctx.fill();
      text(s, (a + b) / 2, y - 11, 13, "#c8402c");
    };
    if (L.key === "a") { bracket(0, 0, "trial state"); bracket(1, 2, "Trotter step"); }
    else { bracket(0, 0, "|111⟩"); bracket(1, G.cols.length - 1, `Trotter step ×${G.n} · T = ${G.T}`); }
  }
  function drawGate(t) {
    const x = xPos(MEASURE_X);
    if (!gateAlive()) { text("gate destroyed", x, BOT + 34, 12, "#8a2a1a"); return; }
    ctx.save();
    if (G.gateFlash > 0) { ctx.shadowColor = "#fff3a0"; ctx.shadowBlur = 50 * G.gateFlash; }
    rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.fillStyle = "#9a6a35"; ctx.fill();
    ctx.restore();
    ctx.lineWidth = 3; ctx.strokeStyle = "#5b3a1a"; rr(x - 12, TOP - 2, 24, DATA_H + 4, 6); ctx.stroke();
    for (const q of lanes) {
      const y = laneY(q);
      circle(x, y, 18); ctx.fillStyle = "#fff8dc"; ctx.fill(); ctx.strokeStyle = "#5b3a1a"; ctx.lineWidth = 2.5; ctx.stroke();
      const ang = -Math.PI / 2 + (G.gateFlash > 0 ? Math.sin(t * 40 + q) * 0.9 : 0.5 * Math.sin(t + q));
      ctx.beginPath(); ctx.moveTo(x, y + 5); ctx.lineTo(x + Math.cos(ang) * 12, y + 5 + Math.sin(ang) * 12);
      ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2.5; ctx.stroke();
    }
    rr(x - 60, BOT + 30, 120, 20, 6); ctx.fillStyle = "#8a5a2b"; ctx.fill();
    text(`Measure · ${G.gateHP} HP`, x, BOT + 40, 12, "#fff4cf");
  }
  function drawRegs() {
    for (const r of G.regs) {
      const X = xPos(r.x);
      for (const q of lanes) {
        const y = laneY(q) - 8;
        ctx.save();
        if (r.measured) {
          circle(X, y, r.flash > 0 ? 10 : 8); ctx.fillStyle = r.basis === "Z" ? "#8be04e" : "#7fd8f0"; ctx.fill();
          ctx.lineWidth = 2; ctx.strokeStyle = "#2f5d1f"; ctx.stroke(); ctx.restore();
          text(sym(r, q), X, y - 17, 13, "#1d3d10");
        } else {
          ctx.shadowColor = "#b388ff"; ctx.shadowBlur = 12;
          circle(X, y, 8); ctx.fillStyle = "#d8c2ff"; ctx.fill(); ctx.lineWidth = 2; ctx.strokeStyle = "#5b3a8c"; ctx.stroke(); ctx.restore();
        }
      }
      if (r.measured) {
        ctx.save(); ctx.setLineDash([4, 4]); ctx.strokeStyle = r.basis === "Z" ? "rgba(47,93,31,.6)" : "rgba(31,111,154,.6)"; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.moveTo(X, laneY(0) - 8); ctx.lineTo(X, laneY(2) - 8); ctx.stroke(); ctx.restore();
        rr(X - 36, laneY(2) + 14, 72, 16, 4); ctx.fillStyle = "rgba(40,25,10,.75)"; ctx.fill();
        text(`${r.basis} avg ${fmtE(r.e)}`, X, laneY(2) + 22, 11, "#fff4cf");
      } else if (G.cols.length) {
        ctx.save(); ctx.setLineDash([5, 4]); ctx.strokeStyle = "rgba(120,70,200,.6)"; ctx.lineWidth = 2;
        ctx.beginPath(); ctx.moveTo(X, laneY(0) - 8); ctx.lineTo(X, laneY(2) - 8); ctx.stroke(); ctx.restore();
      }
    }
  }
  function drawZombieBody(x, y, t, o) {
    const sw = o.still ? 0 : Math.sin(t * 2.6 + o.phase);
    ctx.save(); ctx.translate(x, y); ctx.scale(o.scale, o.scale);
    const coat = o.flash ? "#ffffff" : o.coat, skin = o.flash ? "#ffffff" : o.skin;
    ctx.lineCap = "round"; ctx.strokeStyle = "#3b3f5c"; ctx.lineWidth = 7;
    ctx.beginPath(); ctx.moveTo(-6, 0); ctx.lineTo(-9 + sw * 5, 28); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(6, 0); ctx.lineTo(8 - sw * 5, 28); ctx.stroke();
    rr(-15, -44, 30, 48, 8); ctx.fillStyle = coat; ctx.fill();
    ctx.strokeStyle = skin; ctx.lineWidth = 7;
    ctx.beginPath(); ctx.moveTo(-10, -34); ctx.lineTo(-36, -34 + sw * 3); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(-6, -26); ctx.lineTo(-32, -25 - sw * 3); ctx.stroke();
    circle(-2, -60, 17); ctx.fillStyle = skin; ctx.fill();
    circle(-9, -64, 5.5); ctx.fillStyle = "#fff"; ctx.fill(); circle(4, -63, 4.5); ctx.fill();
    ctx.fillStyle = "#222"; circle(-10, -63, 2.3); ctx.fill(); circle(3, -62, 2.1); ctx.fill();
    ctx.strokeStyle = "#3a2a1a"; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(-12, -51); ctx.lineTo(2, -52); ctx.stroke();
    ctx.fillStyle = "#ffd23f"; ctx.beginPath();
    ctx.moveTo(-16, -74); ctx.lineTo(-12, -92); ctx.lineTo(-4, -80); ctx.lineTo(2, -94); ctx.lineTo(8, -80); ctx.lineTo(14, -90); ctx.lineTo(12, -74);
    ctx.closePath(); ctx.fill(); ctx.strokeStyle = "#9a5b10"; ctx.stroke();
    ctx.restore();
  }
  function drawBoss(t) {
    const B = G.boss;
    if (!B) return;
    const x = xPos(B.x);
    ctx.save(); ctx.strokeStyle = "rgba(224,138,30,.85)"; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.ellipse(x - 4, TOP + DATA_H / 2, 46, DATA_H / 2 - 6, 0, 0, 2 * PI); ctx.stroke(); ctx.restore();
    for (const q of lanes) drawZombieBody(x, laneY(q) + 30, t, { scale: 0.85, phase: q * 1.3, flash: B.flash > 0, still: G.phase === "plan", coat: "#2e3a5e", skin: "#9fc0e8" });
    // Card: the Hamiltonian, then a gauge of the live ⟨H⟩ ± 3σ against the shield.
    const L = lvl(), s = `Hamiltonian boss · ${L.hText}`, P = streamProof();
    const line = P.est
      ? `live ⟨H⟩ ≈ ${fmtE(P.e)} ± ${P.se.toFixed(2)} (${P.n * L.groups.length} shots) · shield ${fmtE(L.bar)}`
        + (P.done ? "" : P.e < L.bar ? (P.n < MIN_SHOTS ? " · need more shots" : " · not sure yet") : "")
      : `live ⟨H⟩: waiting for shots · shield ${fmtE(L.bar)}`;
    ctx.font = `13px ${FONT}`;
    const w1 = ctx.measureText(s).width;
    ctx.font = `12px ${FONT}`;
    const bw = Math.max(260, w1, ctx.measureText(line).width) + 24, bx = Math.min(W - bw / 2 - 6, Math.max(bw / 2 + 6, x));
    rr(bx - bw / 2, 3, bw, 60, 7); ctx.fillStyle = "rgba(30,25,60,.92)"; ctx.fill();
    text(s, bx, 14, 13, "#efe6ff");
    const lo = G.spec[0] - 0.2, hi = G.spec[DIM - 1] + 0.2, gx0 = bx - bw / 2 + 12, gx1 = bx + bw / 2 - 12;
    const gx = e => gx0 + ((Math.max(lo, Math.min(hi, e)) - lo) / (hi - lo)) * (gx1 - gx0);
    rr(gx0, 27, gx1 - gx0, 8, 4); ctx.fillStyle = "#3a3358"; ctx.fill();
    rr(gx0, 27, gx(L.bar) - gx0, 8, 4); ctx.fillStyle = "rgba(182,242,122,.35)"; ctx.fill();
    ctx.fillStyle = "#ff6b5e"; ctx.fillRect(gx(L.bar) - 1.5, 23, 3, 16);
    if (P.est) {
      const col = P.done ? "#b6f27a" : P.e < L.bar ? "#ffd23f" : "#c9a7ff";
      ctx.strokeStyle = col; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(gx(P.e - 3 * P.se), 31); ctx.lineTo(gx(P.e + 3 * P.se), 31); ctx.stroke();
      circle(gx(P.e), 31, 5); ctx.fillStyle = col; ctx.fill();
    }
    text(line, bx, 50, 12, "#efe6ff");
  }
  function drawFx() {
    for (const s of G.splats) {
      ctx.save(); ctx.globalAlpha = Math.max(0, Math.min(1, s.life * 2.5));
      circle(s.x, s.y, s.r); ctx.strokeStyle = s.color; ctx.lineWidth = 3; ctx.stroke(); ctx.restore();
    }
    for (const f of G.floats) {
      ctx.save(); ctx.globalAlpha = Math.min(1, f.life * 2); text(f.s, f.x, f.y, 15, f.good ? "#1d6a0f" : "#5b2a8c"); ctx.restore();
    }
  }
  // ---------- New-boss callouts ----------
  // The first time each Hamiltonian boss shows up, a parchment tag with a red bracket (the same
  // style as the circuit tags) says what it is and what beats it. It stays up while you plan and
  // for CALLOUT_SECONDS of battle time.
  const CALLOUT_SECONDS = 7;
  const CALLOUTS = {
    a: { title: "Hamiltonian boss", lines: [
      "Its shield is an energy rule: H = X₀X₁ + X₁X₂.",
      "Send it low-energy states. Your shots must prove",
      "the energy ⟨H⟩ is below −1.8 to break the shield."] },
    b: { title: "Field boss", lines: [
      "A new part ½(Z₀ + Z₁ + Z₂) pulls against the X·X part.",
      "Slide down to its ground state with Trotter steps,",
      "and prove ⟨H⟩ is below −2.2."] }
  };
  const zombieKind = () => null;
  const bossKind = () => lvl().key;
  const zombieAnchor = () => null;
  const bossSpan = () => ({ x: xPos(G.boss.x) - 62, y1: TOP + 8, y2: TOP + DATA_H - 8 });
  const seenKinds = new Set();
  let callouts = [], calloutClock = 0;
  function meetKind(kind, target, boss) {
    if (!kind || seenKinds.has(kind)) return;
    seenKinds.add(kind);
    callouts.push({ kind, target, boss, left: CALLOUT_SECONDS });
  }
  function tagSize(C) {
    ctx.font = `14px ${FONT}`;
    let w = ctx.measureText(`NEW · ${C.title}`).width;
    ctx.font = `12px ${FONT}`;
    for (const s of C.lines) w = Math.max(w, ctx.measureText(s).width);
    return { w: w + 20, h: 24 + 16 * C.lines.length };
  }
  function drawTag(C, x, y, m) {
    rr(x, y, m.w, m.h, 6); ctx.fillStyle = "rgba(255,248,220,.96)"; ctx.fill();
    ctx.strokeStyle = "rgba(200,64,44,.55)"; ctx.lineWidth = 1.5; ctx.stroke();
    text(`NEW · ${C.title}`, x + 10, y + 13, 14, "#c8402c", "left");
    C.lines.forEach((s, i) => text(s, x + 10, y + 30 + 16 * i, 12, "#5b3a1a", "left"));
  }
  // Tag above the zombie (below it when there is no room), with a bracket over its head.
  function zombieCallout(c) {
    if (!c.gone) c.a = zombieAnchor(c.target);
    const C = CALLOUTS[c.kind], a = c.a, m = tagSize(C);
    const above = a.top - 22 - m.h >= 4, d = above ? 6 : -6;
    const by = above ? a.top - 8 : a.bottom + 8, ty = above ? by - 14 - m.h : by + 14;
    const tx = Math.max(4, Math.min(W - m.w - 4, a.x - m.w / 2));
    ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(a.x - 26, by + d); ctx.quadraticCurveTo(a.x - 26, by, a.x - 18, by);
    ctx.lineTo(a.x + 18, by); ctx.quadraticCurveTo(a.x + 26, by, a.x + 26, by + d); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(a.x, by); ctx.lineTo(Math.max(tx + 12, Math.min(tx + m.w - 12, a.x)), above ? ty + m.h : ty); ctx.stroke();
    drawTag(C, tx, ty, m);
  }
  // Tag to the left of the boss, with a bracket along its whole height.
  function bossCallout(c) {
    if (!c.gone) c.s = bossSpan();
    const C = CALLOUTS[c.kind], s = c.s, m = tagSize(C), ym = (s.y1 + s.y2) / 2;
    const tx = Math.max(4, s.x - 16 - m.w), ty = Math.max(4, ym - m.h / 2);
    ctx.strokeStyle = "#c8402c"; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(s.x + 6, s.y1); ctx.quadraticCurveTo(s.x, s.y1, s.x, s.y1 + 8);
    ctx.lineTo(s.x, s.y2 - 8); ctx.quadraticCurveTo(s.x, s.y2, s.x + 6, s.y2); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(s.x, ym); ctx.lineTo(tx + m.w, ym); ctx.stroke();
    drawTag(C, tx, ty, m);
  }
  function drawCallouts(t) {
    const dt = Math.min(0.05, Math.max(0, t - (calloutClock || t)));
    calloutClock = t;
    for (const z of G.zombies || []) if (z.alive && !(z.wait > 0)) meetKind(zombieKind(z), z, false);
    if (G.boss) meetKind(bossKind(), G.boss, true);
    // A tag whose zombie has just fallen stays where it was for up to 3 more seconds.
    for (const c of callouts) {
      if (!c.gone && !(c.boss ? G.boss === c.target : c.target.alive && G.zombies.includes(c.target))) { c.gone = true; c.left = Math.min(c.left, 3); }
    }
    callouts = callouts.filter(c => c.left > 0);
    for (const c of callouts) {
      if ((G.phase === "battle" && !G.paused) || c.gone) c.left -= dt;
      ctx.save(); ctx.globalAlpha = Math.max(0, Math.min(1, c.left / 0.6));
      if (c.boss) bossCallout(c); else zombieCallout(c);
      ctx.restore();
    }
  }

  function draw(t) {
    ctx.clearRect(0, 0, W, H);
    drawLawn(); drawHouse(); drawCircuit(t); drawGate(t); drawBoss(t); drawRegs(); drawFx(); drawCallouts(t);
    if (G.phase === "plan") {
      const a = xPos(6.0), b = xPos(12.4);
      rr(a, TOP + 30, b - a, 34, 8); ctx.fillStyle = "rgba(40,25,10,.8)"; ctx.fill();
      text(`Level ${lvl().key}: ${lvl().name}. Set your circuit, then press Ready.`, (a + b) / 2, TOP + 47, 14, "#fff4cf");
    }
    if (G.paused || G.phase === "over") {
      ctx.fillStyle = "rgba(20,40,10,.55)"; ctx.fillRect(0, 0, W, H);
      text(G.phase === "over" ? "Base lost" : "Paused", W / 2, H / 2 - 10, 40, "#fff4cf");
      text(G.phase === "over" ? "Restart to play again" : "Resume when ready", W / 2, H / 2 + 26, 16, "#fff4cf");
    }
  }

  // ---------- Events ----------
  document.querySelectorAll("#levelSel button").forEach((b, i) => b.addEventListener("click", () => {
    if (G.phase === "battle" && !confirm("Leave the level in progress?")) return;
    loadLevel(i, false);
  }));
  $("zeroP").addEventListener("click", () => { Object.keys(G.p).forEach(k => { G.p[k] = 0; }); circuitChanged(); say("All angles reset to 0."); });
  $("preview").addEventListener("change", e => { G.preview = e.target.checked; drawSpectrum(); drawLandscape(); updateHUD(); });
  $("ready").addEventListener("click", startWave);
  $("duel").addEventListener("click", duel);
  $("pause").addEventListener("click", () => { if (G.phase === "battle") { G.paused = !G.paused; updateHUD(); } });
  $("restart").addEventListener("click", () => loadLevel(0, false));
  canvas.addEventListener("click", ev => {
    const rc = canvas.getBoundingClientRect();
    const x = (ev.clientX - rc.left) * (W / rc.width), y = (ev.clientY - rc.top) * (H / rc.height);
    const col = G.cols.find(c => Math.abs(xPos(c.u) - x) < 26);
    if (!col) return;
    const tl = col.tiles.find(t => y >= laneY(Math.min(...t.qs)) - 30 && y <= laneY(Math.max(...t.qs)) + 30);
    if (!tl) return;
    if (tl.locked || !tl.param) { say(lvl().key === "b" ? "These angles follow from T, n and the order. Use the controls below." : "This gate is locked in place."); return; }
    G.sel = tl.param; syncParams(); drawLandscape();
    say(`Selected ${tl.param}. Drag its slider below${G.preview ? "; the landscape shows E against it" : ""}.`);
  });

  // ---------- Tutorial storybook ----------
  // Pages shown before each level: a kid-friendly explanation, a "Try this" box and a collapsible
  // Physicist's corner. gt() draws a gate the way it looks on the lawn.
  const gt = (s, k) => `<span class="gt gt-${k}">${s}</span>`;
  const TUT = {
    intro: [{
      title: "Lowest energy wins",
      kid: `<p>A ball always rolls down to the bottom of a valley. Nature loves the <b>lowest energy</b>, and quantum things do too. Their lowest-energy state is called the <b>ground state</b>.</p>
        <p>The <b>Hamiltonian boss</b> carries an energy shield. Its rule for energy, the <b>Hamiltonian H</b>, is written on its card. Build a state whose energy is <b>below the shield</b> and the shield breaks!</p>`,
      phys: `<p>This is the idea behind the variational quantum eigensolver (VQE): prepare a trial state |ψ(θ)⟩ and estimate ⟨ψ|H|ψ⟩ from measurements. By the variational principle ⟨H⟩ ≥ E₀, so a lower estimate is always a better approximation to the ground state.</p>`
    }, {
      title: "Measuring energy with shots",
      kid: `<p>You can't read energy straight off a qubit. Instead every register fires a <b>volley of shots</b> (6 in level a, 200 in level b), and the game <b>averages</b> them, like rolling a die many times to find its average. More shots give a sharper answer.</p>
        <p>The gauge above the boss shows the average (the dot), how unsure it still is (the whisker) and the shield (the red line). When the <b>whole whisker</b> is past the red line, the boss falls. Impatient? Press <b>Energy duel</b> to fire 8000 shots at once.</p>`,
      try: ["Watch the dot settle down as more shots arrive."],
      phys: `<p>H is split into groups that can be measured together: X-basis shots for the X·X terms and Z-basis shots for the Z terms. ⟨H⟩ is estimated as the sum of the group means; the whisker is ±3 standard errors, and the stream defeats the boss once estimate + 3σ < shield (about 99.7% confidence).</p>`
    }],
    levels: [
      [{
        title: "The boss's energy rule",
        kid: `<p>The rule is <b class="bk">H = X₀X₁ + X₁X₂</b>. Each X·X part looks at two neighbouring qubits along the sideways (x) direction. Pointing the <b>same way</b> (+ +) costs +1. Pointing <b>opposite ways</b> (+ −) gives −1.</p>
          <p>So the lowest energy, <b>−2</b>, comes when neighbours take turns: <b>+ − +</b> or <b>− + −</b>. The shield is at −1.8.</p>`,
        phys: `<p>Both terms are diagonal in the X basis, so E = x₀x₁ + x₁x₂ with xᵢ = ±1. The ground space at −2 is two-fold degenerate: |+−+⟩ and |−+−⟩. (This is X₁X₂ + X₂X₃ in 1-based labels.)</p>`
      }, {
        title: "Aiming with RY",
        kid: `<p>${gt("RY(θ)", "r")} aims each qubit. θ = π/2 (1.571) points it to <b>+</b>. θ = 3π/2 (4.712) points it to <b>−</b>. θ = 0 points straight up at |0⟩, exactly halfway between + and − along x.</p>
          <p>Drag the sliders and watch the live energy fall.</p>`,
        try: ["θ₀ = 4.712, θ₁ = 1.571, θ₂ = 4.712 makes − + − with energy −2.",
          "Click a gate on the lawn to pick its slider.",
          "Tick <b>Simulator preview</b> to see the energy curve for the selected angle."],
        phys: `<p>⟨X⟩ for RY(θ)|0⟩ is sin θ, so before the RXX step E(θ) = sin θ₀ sin θ₁ + sin θ₁ sin θ₂, minimised at θ₁ = π/2 with θ₀ = θ₂ = 3π/2, or the mirror choice.</p>`
      }, {
        title: "The time step that changes nothing",
        kid: `<p>The two ${gt("RXX", "c")} gates are one <b>time step</b> of the boss's own rule: a <b>Trotter step</b>. Spinning a ball that sits at the bottom of a bowl doesn't lift it out, so the energy stays the same.</p>
          <p>Try moving φ: the energy never changes!</p>`,
        phys: `<p>Because the two terms commute, RXX(φ)₀₁ · RXX(φ)₁₂ = e^(−iφH/2) exactly. Since [H, e^(−iHt)] = 0, ⟨H⟩ is conserved under H's own time evolution. Lowering the energy needs an operation that does not commute with H, here the RY rotations.</p>`
      }],
      [{
        title: "A trickier rule",
        kid: `<p>The boss adds a new part: <b class="bk">½(Z₀ + Z₁ + Z₂)</b>. This part wants every qubit pointing <b>down</b>, to |1⟩. The X·X part wants them pointing sideways. The two parts <b>pull against each other</b>, so the best state is a quantum mixture of both and hard to guess by hand.</p>`,
        phys: `<p>[X₀X₁, Z₀] ≠ 0, so the terms no longer commute. The ground energy is −2.403 and the first excited level −2.209. The ground state is entangled and lies in the odd Z-parity sector, the same sector as |111⟩.</p>`
      }, {
        title: "Change slowly: the adiabatic trick",
        kid: `<p>Start with the easy part alone: all qubits down, <b class="bk">|111⟩</b>, the ground state of the Z part. Then <b>slowly</b> turn up the X·X part.</p>
          <p>It's like carrying a sleeping cat to a new bed. Move slowly and it stays asleep: the qubits stay in the lowest-energy state. Rush, and the cat wakes up. The total time you take is <b>T</b>.</p>`,
        phys: `<p>H(s) = s(X₀X₁ + X₁X₂) + (1 − s/2)ΣZ, with s going from 0 to 1 over time T. By the adiabatic theorem the state follows the instantaneous ground state when T is long compared with the inverse square of the smallest energy gap along the path.</p>`
      }, {
        title: "Trotter–Suzuki: slicing time",
        kid: `<p>A quantum computer can't apply both parts at the same moment, so it <b>takes turns</b>: a little ${gt("RZ", "z")}, a little ${gt("RXX", "c")}, a little RZ, a little RXX... like walking left foot, right foot. Each turn is one <b>Trotter step</b>, and <b>n</b> is how many steps you take.</p>
          <p>More, smaller steps follow the true path more closely. The <b>2nd-order (Suzuki)</b> step goes half RZ, full RXX, half RZ. It's balanced, and much more accurate.</p>`,
        try: ["Try 2nd order, n = 12, T = 5: the boss falls in about 5 registers.",
          "Tick <b>Simulator preview</b>: the dashed curve is perfect evolution with no slicing, the solid curve is yours.",
          "T too small: the cat wakes up. Too few steps: slicing errors. The shield is −2.20."],
        phys: `<p>Lie–Trotter: e^(−iAΔt) e^(−iBΔt) = e^(−i(A+B)Δt) + O(Δt²) per step. The symmetric Suzuki step e^(−iAΔt/2) e^(−iBΔt) e^(−iAΔt/2) has O(Δt³) error per step, so far fewer steps reach the same accuracy.</p>`
      }]
    ]
  };
  const tut = { pages: [], i: 0, seen: new Set(), resume: false };
  function openTutorial(withIntro) {
    const tag = `Level ${lvl().key}`;
    tut.pages = [...(withIntro ? TUT.intro.map(p => ({ ...p, tag: "Welcome" })) : []),
      ...(TUT.levels[G.level] || []).map(p => ({ ...p, tag }))];
    if (!tut.pages.length) return;
    tut.i = 0;
    tut.resume = G.phase === "battle" && !G.paused;
    if (tut.resume) G.paused = true;
    renderTutorial();
    $("tut").classList.remove("hidden");
    $("tutNext").focus();
  }
  function closeTutorial() {
    $("tut").classList.add("hidden");
    if (tut.resume) { G.paused = false; tut.resume = false; }
    updateHUD();
  }
  function renderTutorial() {
    const p = tut.pages[tut.i], n = tut.pages.length;
    $("tutStep").textContent = `${p.tag} · page ${tut.i + 1} of ${n}`;
    $("tutTitle").textContent = p.title;
    $("tutBody").innerHTML = p.kid + (p.try ? `<div class="tut-try"><b>Try this</b><ul>${p.try.map(t => `<li>${t}</li>`).join("")}</ul></div>` : "");
    $("tutPhysWrap").classList.toggle("hidden", !p.phys);
    $("tutPhys").innerHTML = p.phys || "";
    $("tutPhysWrap").open = false;
    $("tutBack").disabled = tut.i === 0;
    $("tutNext").textContent = tut.i === n - 1 ? "Let's play!" : "Next";
    $("tutDots").innerHTML = tut.pages.map((_, k) => `<i class="${k === tut.i ? "on" : ""}"></i>`).join("");
  }
  // Each level's pages open once per session, the first time it loads (with the welcome pages before level a).
  function autoTutorial() {
    if (tut.seen.has(G.level)) return;
    const first = tut.seen.size === 0;
    tut.seen.add(G.level);
    openTutorial(first);
  }
  $("tutNext").addEventListener("click", () => { if (tut.i < tut.pages.length - 1) { tut.i++; renderTutorial(); } else closeTutorial(); });
  $("tutBack").addEventListener("click", () => { if (tut.i > 0) { tut.i--; renderTutorial(); } });
  $("tutSkip").addEventListener("click", closeTutorial);
  $("tutOpen").addEventListener("click", () => openTutorial(true));
  document.addEventListener("keydown", e => {
    if ($("tut").classList.contains("hidden")) return;
    if (e.key === "Escape") closeTutorial();
    else if (e.key === "ArrowRight") $("tutNext").click();
    else if (e.key === "ArrowLeft") $("tutBack").click();
  });

  // ---------- Loop ----------
  let last = 0;
  function frame(ms) {
    const t = ms / 1000, dt = Math.min(t - (last || t), 0.05);
    last = t;
    update(dt);
    draw(t);
    G.monT += dt;
    if (G.monT > 0.25) { G.monT = 0; updateHUD(); drawSpectrum(); }
    requestAnimationFrame(frame);
  }

  window.QPVZ3 = { G, LEVELS, update, loadLevel, startWave, duel, circuitChanged, exactE, spectrum, adiabatic, energy, basisState, HFIELD, XX };
  fitCanvas();
  window.addEventListener("resize", () => { fitCanvas(); drawSpectrum(); drawLandscape(); });
  loadLevel(0, false);
  requestAnimationFrame(frame);
})();
