// Quantum Battlefield launcher: home (three floating adventure tiles) and, at #adv-N, that
// adventure's level map. Levels open in "<folder>/play.html?level=i".
(function () {
  const P = window.QBProgress;
  const $ = id => document.getElementById(id);
  const SVGNS = "http://www.w3.org/2000/svg";

  // Map positions are in the board's own units: wide board 160 x 100, tall (phone) board 100 x 130.
  const ADVENTURES = {
    1: {
      name: "Qubit basics", folder: "1st%20adventure", color: "var(--a1)",
      wide: [[30, 70], [80, 32], [130, 64]],
      tall: [[32, 104], [68, 66], [34, 26]],
      levels: [
        { key: "A", name: "First qubit", desc: "Flip and split one qubit with X, H and Z, then measure it in the Z or X basis to hit the classical zombies." },
        { key: "B", name: "Rotations", desc: "Turn the RX, RY and RZ dials until your bullet copies a quantum zombie's state exactly. Only a perfect match makes it vanish." },
        { key: "C", name: "Entangled boss", desc: "Two qubits and a CNOT. The Bell-state boss only falls to an entangled pair; no ordinary pair can prove itself." }
      ]
    },
    2: {
      name: "Grover's search", folder: "2nd%20adventure", color: "var(--a2)",
      wide: [[24, 72], [60, 30], [100, 70], [138, 32]],
      tall: [[28, 112], [72, 84], [30, 52], [70, 20]],
      levels: [
        { key: "a", name: "Mark the boss", desc: "Build the oracle. Phase kickback flips the boss's hiding pattern negative. Beat a 50% hit rate." },
        { key: "b", name: "Build the diffuser", desc: "Add the diffuser so interference grows the marked pattern and shrinks the rest. Beat 70%." },
        { key: "c", name: "Superposed boss", desc: "The boss hides in two patterns at once and standard Grover freezes at 50%. Change the marking phase to break it." },
        { key: "d", name: "Three lanes", desc: "Eight patterns on three qubits. Repeat Grover's block just the right number of times to beat 90%." }
      ]
    },
    3: {
      name: "Ground-state hunters", folder: "3rd%20adventure", color: "var(--a3)",
      wide: [[42, 66], [118, 36]],
      tall: [[34, 96], [66, 34]],
      levels: [
        { key: "a", name: "Commuting ladder", desc: "Explore the rotation angles to drive the energy of H = X₀X₁ + X₁X₂ below the boss's −1.8 shield." },
        { key: "b", name: "Trotter–Suzuki", desc: "Switch the couplings on slowly with repeated Trotter steps and slide down to the ground state, below −2.2." }
      ]
    }
  };

  const STAR = '<svg viewBox="0 0 24 24"><path class="CLS" d="M12 2.6l2.8 5.9 6.4.8-4.7 4.4 1.2 6.4L12 17l-5.7 3.1 1.2-6.4L2.8 9.3l6.4-.8z"/></svg>';
  const starSvg = on => STAR.replace("CLS", on ? "star-on" : "star-off");
  const LOCK_CLOSED = '<rect x="9" y="21" width="30" height="22" rx="4" fill="#ffd23f" stroke="#5b3a1a" stroke-width="3.5"/>' +
    '<path d="M15 21v-6a9 9 0 0 1 18 0v6" fill="none" stroke="#5b3a1a" stroke-width="3.5"/>' +
    '<circle cx="24" cy="30" r="3" fill="#5b3a1a"/><path d="M24 31v6" stroke="#5b3a1a" stroke-width="3" stroke-linecap="round"/>';
  const LOCK_OPEN = '<rect x="9" y="21" width="30" height="22" rx="4" fill="#ffd23f" stroke="#5b3a1a" stroke-width="3.5"/>' +
    '<path d="M15 21v-6a9 9 0 0 1 17.6-2.6" fill="none" stroke="#5b3a1a" stroke-width="3.5" stroke-linecap="round"/>' +
    '<circle cx="24" cy="30" r="3" fill="#5b3a1a"/><path d="M24 31v6" stroke="#5b3a1a" stroke-width="3" stroke-linecap="round"/>';
  const lockIcon = open => `<svg viewBox="0 0 48 48">${open ? LOCK_OPEN : LOCK_CLOSED}</svg>`;

  const phone = window.matchMedia("(max-width: 640px)");
  let advId = null, sel = 0;

  // ---------- Home ----------
  function renderHome() {
    document.querySelectorAll(".tile-stars").forEach(el => {
      const a = +el.dataset.adv, n = ADVENTURES[a].levels.length;
      let got = 0;
      for (let i = 0; i < n; i++) got += P.stars(a, i);
      el.innerHTML = `<b>★</b> ${got} / ${3 * n}`;
    });
  }

  // ---------- Map ----------
  // A wobbly island outline around (cx, cy), smoothed through edge midpoints.
  function blobPath(cx, cy, r, seed) {
    const k = 9, pts = [];
    for (let i = 0; i < k; i++) {
      const a = (i / k) * Math.PI * 2;
      const rr = r * (1 + 0.16 * Math.sin(3 * a + seed) + 0.1 * Math.cos(5 * a + 2 * seed));
      pts.push([cx + rr * Math.cos(a), cy + rr * Math.sin(a)]);
    }
    const mid = (p, q) => [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2];
    let d = `M${mid(pts[k - 1], pts[0]).map(v => v.toFixed(2)).join(" ")}`;
    for (let i = 0; i < k; i++) {
      const p = pts[i], m = mid(p, pts[(i + 1) % k]);
      d += ` Q${p[0].toFixed(2)} ${p[1].toFixed(2)} ${m[0].toFixed(2)} ${m[1].toFixed(2)}`;
    }
    return d + "Z";
  }

  function svgEl(tag, attrs) {
    const el = document.createElementNS(SVGNS, tag);
    for (const k in attrs) el.setAttribute(k, attrs[k]);
    return el;
  }

  function renderMap() {
    const A = ADVENTURES[advId], tall = phone.matches;
    const W = tall ? 100 : 160, H = tall ? 130 : 100;
    const pos = tall ? A.tall : A.wide;
    const n = A.levels.length;

    $("mapNum").textContent = `Adventure ${advId}`;
    $("mapTitle").textContent = A.name;
    document.title = `Quantum Battlefield · ${A.name}`;

    const svg = $("boardSvg");
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.setAttribute("preserveAspectRatio", "none");
    svg.replaceChildren();
    pos.forEach(([x, y], i) => svg.appendChild(svgEl("path", { class: "blob", d: blobPath(x, y, tall ? 20 : 19, i * 1.7 + advId) })));
    for (let i = 1; i < n; i++) {
      const [x1, y1] = pos[i - 1], [x2, y2] = pos[i];
      const open = P.unlocked(advId, i);
      if (open) svg.appendChild(svgEl("line", { class: "trail-edge", x1, y1, x2, y2 }));
      svg.appendChild(svgEl("line", { class: "trail" + (open ? "" : " locked"), x1, y1, x2, y2 }));
    }

    const box = $("nodes");
    box.replaceChildren();
    A.levels.forEach((L, i) => {
      const open = P.unlocked(advId, i), s = P.stars(advId, i);
      const b = document.createElement("button");
      b.type = "button";
      b.className = "node" + (open ? "" : " locked") + (i === sel ? " sel" : "");
      b.style.left = (pos[i][0] / W * 100) + "%";
      b.style.top = (pos[i][1] / H * 100) + "%";
      b.setAttribute("role", "listitem");
      b.setAttribute("aria-pressed", i === sel ? "true" : "false");
      b.setAttribute("aria-label", `Level ${L.key}, ${L.name}: ${open ? `${s} of 3 stars` : "locked"}`);
      b.innerHTML = `<span class="disc">${open ? L.key : lockIcon(false)}</span>` +
        (open ? `<span class="stars">${[0, 1, 2].map(k => starSvg(k < s)).join("")}</span>` : "");
      b.addEventListener("click", () => { sel = i; renderMap(); });
      b.addEventListener("dblclick", () => { if (open) location.href = playUrl(i); });
      box.appendChild(b);
    });

    // Star chart
    $("bars").innerHTML = A.levels.map((L, i) => {
      const s = P.stars(advId, i);
      return `<div class="bar${s ? "" : " zero"}"><span style="height:${s ? (s / 3) * 100 : 4}%"></span></div>`;
    }).join("");
    let labels = $("bars").nextElementSibling;
    if (!labels || !labels.classList.contains("bar-labels")) {
      labels = document.createElement("div");
      labels.className = "bar-labels";
      $("bars").after(labels);
    }
    labels.innerHTML = A.levels.map(L => `<span>${L.key}</span>`).join("");

    // Door: the selected level
    const L = A.levels[sel], open = P.unlocked(advId, sel), s = P.stars(advId, sel);
    $("padlock").innerHTML = open ? LOCK_OPEN : LOCK_CLOSED;
    $("doorKey").textContent = `Level ${L.key}` + (open && s ? ` · ${"★".repeat(s)}` : "");
    $("doorName").textContent = L.name;
    $("doorDesc").textContent = open ? L.desc : `Locked. Clear level ${A.levels[sel - 1].key} to unlock it.`;
    const play = $("play");
    play.href = open ? playUrl(sel) : "#";
    play.textContent = open ? (s ? "Play again" : "Play") : "Locked";
    play.setAttribute("aria-disabled", open ? "false" : "true");
    play.tabIndex = open ? 0 : -1;

    $("unlockAll").checked = P.all;
  }

  const playUrl = i => `${ADVENTURES[advId].folder}/play.html?level=${i}`;

  // Pick the level to highlight: the first open level without stars, else the last open one.
  function defaultSel(a) {
    const n = ADVENTURES[a].levels.length;
    let last = 0;
    for (let i = 0; i < n; i++) {
      if (!P.unlocked(a, i)) break;
      last = i;
      if (!P.stars(a, i)) return i;
    }
    return last;
  }

  // ---------- Routing ----------
  function route(focus) {
    P.refresh();
    const m = /^#adv-([123])$/.exec(location.hash);
    if (m) {
      const a = +m[1];
      if (a !== advId) { advId = a; sel = defaultSel(a); }
      $("home").hidden = true;
      $("map").hidden = false;
      renderMap();
      if (focus) $("mapTitle").focus();
    } else {
      advId = null;
      document.title = "Quantum Battlefield";
      $("map").hidden = true;
      $("home").hidden = false;
      renderHome();
    }
    if (focus) window.scrollTo(0, 0);
  }

  window.addEventListener("hashchange", () => route(true));
  window.addEventListener("pageshow", () => route(false));
  phone.addEventListener("change", () => { if (advId) renderMap(); });
  $("unlockAll").addEventListener("change", e => { P.setAll(e.target.checked); renderMap(); });
  $("reset").addEventListener("click", () => {
    if (!confirm("Clear all stars and lock the levels again?")) return;
    P.reset();
    sel = 0;
    renderMap();
  });
  route(false);
})();
