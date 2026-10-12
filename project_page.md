# Quantum Battlefield: a quantum twist on Plants vs. Zombies

**Your bullets are qubits. Your plants are quantum gates. The zombies obey quantum mechanics.**

> 🎬 **Demo video:** *(link coming soon)*

---

## The idea

Every shot fired across the lawn is a **live quantum state**. It stays quantum, in superposition, carrying amplitudes and phases, until something measures it and it **collapses** into an ordinary classical pea. Between the launcher and the battlefield sits your quantum circuit: the gates you plant are your **powers to manipulate the state**, and the zombies only fall to the state you actually built.

The game is a quantum circuit you can watch in motion.

---

## What each part of the game means

| In the game | In quantum physics |
|---|---|
| **Bullets** flying down a lane | **Qubits**: live quantum states, until they are measured |
| **Plants** in the slots | **Quantum gates**: your powers to rotate, flip, phase-shift and entangle the state |
| **Number of slots** | **Circuit depth**. We deliberately limit it, so every gate has to count, just as on real hardware where depth costs accuracy |
| **The wooden gate** mid-lawn | **Measurement**. The bullet collapses, and a classical answer is printed on it |
| **Classical zombies** | Fall only to **measured, classical** bullets carrying the right answer |
| **Quantum zombies** | Carry a quantum state of their own and fall only to an **unmeasured bullet in the same state**. Classical bullets pass straight through them, and so do quantum bullets in a different state |
| **The entangled boss** | A Bell state that only an **entangled** pair of bullets can defeat |
| **The Hamiltonian boss** | An energy rule, defeated only by a **low-energy state** |

---

## Watch the quantum wave

As each bullet passes a gate, the game shows its **probability distribution and amplitudes updating live**: you see superposition spread, phases flip sign, and probabilities shift gate by gate. A Bloch sphere shows the qubit's arrow turning, and a dashed link appears between bullets the moment they become entangled.

The circuit stops being boxes on a diagram. It becomes a wave you watch being shaped.

---

## Adventure 1: Qubit basics

Superposition, measurement bases, phase, rotations and entanglement. You learn to **read a quantum zombie's state** (say, √3/2|0⟩ + i/2|1⟩), recreate it with rotation gates, and see that only the matching quantum bullet makes it vanish. The finale is an entangled boss that **no un-entangled pair can beat**, however good the aim: the game turns the entanglement bound into a wall only real entanglement can break.

## Adventure 2: Seeing Grover's search work

The boss hides in a secret bit pattern. You build Grover's machine (oracle, phase kickback, diffuser) and **watch it work**. The amplitude bars show the marked pattern flip negative, then grow while every other pattern shrinks. Repeat the block too often and you watch it **overshoot**. When the boss hides in two patterns at once, standard Grover freezes at exactly 50%, and you crack it by changing the phase you mark with, from a half turn to a quarter turn.

## Adventure 3: Hamiltonian simulation and VQE

The boss is a Hamiltonian, and its shield breaks only below a target energy. You **explore the parameter space** of your circuit's rotation angles, watching the measured energy fall toward the ground state, exactly as in a variational quantum eigensolver.

Then you **simulate the Hamiltonian itself**. Time evolution is built from a short block of gates (a **Trotter step**) **repeated n times**. You see the trade-off directly: more repetitions and the 2nd-order Suzuki step track the true evolution more accurately, while too few repetitions, or rushing the evolution, leave the state stuck at higher energy. A live landscape compares your circuit with perfect evolution, so the **accuracy of the repeated block** becomes visible.

---

## From game to real code

Whatever you plant, the game writes out the **equivalent Qiskit circuit** for you: gates, measurements and, in Adventure 3, the Hamiltonian as a `SparsePauliOp`. The intuition you build on the lawn transfers directly to a real quantum program.

Every outcome comes from the exact quantum state through the Born rule, and we checked the game against textbook results (Grover's success rates, the entanglement bound, ground-state energies). If you win, the quantum mechanics really did the work.

---

## References

**Quantum computing foundations and Grover's search**
- M. A. Nielsen and I. L. Chuang, *Quantum Computation and Quantum Information*, Cambridge University Press (2010).
- NPTEL, quantum computing courses.
- QURECA, ADEQUATE programme.
- QWorld, QBronze, QSilver and QNickel workshops.

**Hamiltonian simulation and Trotter–Suzuki evolution**
- H. F. Trotter, *On the product of semi-groups of operators*, Proceedings of the American Mathematical Society **10**, 545 (1959). *The first-order product formula behind each Trotter step.*
- M. Suzuki, *Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems*, Communications in Mathematical Physics **51**, 183 (1976). *Systematic product-formula approximations.*
- G. Strang, *On the construction and comparison of difference schemes*, SIAM Journal on Numerical Analysis **5**, 506 (1968). *The symmetric half–full–half splitting used as the second-order step.*
- M. Suzuki, *Fractal decomposition of exponential operators with applications to many-body theories and Monte Carlo simulations*, Physics Letters A **146**, 319 (1990). *Symmetric second-order and higher-order Suzuki formulas.*
- S. Lloyd, *Universal quantum simulators*, Science **273**, 1073 (1996). *Simulating local Hamiltonians on a quantum computer with product formulas.*
- N. Wiebe, D. W. Berry, P. Høyer and B. C. Sanders, *Higher order decompositions of ordered operator exponentials*, Journal of Physics A **43**, 065203 (2010). *Product formulas for time-dependent Hamiltonians like our adiabatic path H(s).*
- A. M. Childs, Y. Su, M. C. Tran, N. Wiebe and S. Zhu, *Theory of Trotter error with commutator scaling*, Physical Review X **11**, 011020 (2021). *Why the Trotter error comes from non-commuting terms, and how it shrinks with more steps.*

**Adiabatic state preparation**
- E. Farhi, J. Goldstone, S. Gutmann and M. Sipser, *Quantum computation by adiabatic evolution*, arXiv:quant-ph/0001106 (2000).
- T. Albash and D. A. Lidar, *Adiabatic quantum computation*, Reviews of Modern Physics **90**, 015002 (2018).

---

## About me

I'm **Sriharsha**, a fourth-year BS-MS physics student at **IISER Kolkata**. My research has taken me from material simulations and the theory of spintronics and orbitronics (IISER Kolkata) to quantum simulation of quantum field theories and Hamiltonian simulation (IIT Hyderabad). Along the way I've trained through NPTEL's quantum computing courses, QWorld's QBronze, QSilver and QNickel workshops, Amaravati Quantum Valley, QURECA's ADEQUATE programme, IKAA's quantum error correction course and QGSS 2026.

I built Quantum Battlefield because the ideas that make quantum computing powerful (superposition, interference, entanglement and Hamiltonian simulation) are usually hidden behind equations. I wanted a place where anyone, from a curious kid to a physics student, can *watch* a quantum state being shaped and feel why it behaves the way it does.

---

## Licence

Released under the **MIT License**. Created by **Sriharsha** for the **Qollab Hackathon 2026**.
