// Quantum Battlefield progress: stars per level, kept in this browser's localStorage.
// Level 0 of every adventure is always open; each later level opens once the one before it
// has at least one star, or when "Unlock all" is switched on. Storage may be unavailable
// (private windows, blocked site data), so everything also works from memory.
window.QBProgress = (function () {
  const KEY = "qb.progress.v1";
  function read() {
    try { return JSON.parse(localStorage.getItem(KEY)) || {}; } catch (e) { return {}; }
  }
  function write() {
    try { localStorage.setItem(KEY, JSON.stringify(mem)); } catch (e) { /* memory only */ }
  }
  let mem = read();
  const stars = (adv, i) => (mem.stars && mem.stars[adv + "-" + i]) || 0;
  return {
    stars,
    unlocked: (adv, i) => i === 0 || !!mem.all || stars(adv, i - 1) > 0,
    record(adv, i, s) {
      mem.stars = mem.stars || {};
      if (s > stars(adv, i)) { mem.stars[adv + "-" + i] = s; write(); }
    },
    get all() { return !!mem.all; },
    setAll(on) { mem.all = !!on; write(); },
    reset() { mem = {}; write(); },
    refresh() { mem = read(); },
    // Stars from the base health left when the level is cleared.
    starsFor: baseHP => (baseHP >= 10 ? 3 : baseHP >= 6 ? 2 : 1)
  };
})();
