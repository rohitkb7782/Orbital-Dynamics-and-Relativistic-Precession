# Orbital Dynamics and Relativistic Precession

A computational physics project using numerical simulation to study orbital dynamics. The project compares Euler's Method, fourth-order Runge-Kutta (RK4), and Velocity Verlet for Newtonian orbits. Analytical Keplerian orbits are used as a reference to measure numerical error and convergence. Long-term energy behavior is also compared between the different methods. A relativistic correction is then introduced to study periapsis precession and the stability of circular orbits near the innermost stable circular orbit (ISCO).

![Relativistic orbital dynamics](images/relativistic_orbits.png?raw=true)

**Figure 3.** *Relativistic orbital dynamics for weak-field, intermediate-radius, and near-ISCO orbits. Each panel compares outward and inward radial perturbations.*

## Table of Contents

* [Motivation](#motivation)
* [Mathematical Model](#mathematical-model)

  * [Newtonian Orbital Dynamics](#newtonian-orbital-dynamics)
  * [Keplerian Analytical Solution](#keplerian-analytical-solution)
  * [Relativistic Correction](#relativistic-correction)
  * [Circular Orbits, Stability, and the ISCO](#circular-orbits-stability-and-the-isco)
  * [Relativistic Precession](#relativistic-precession)
* [Numerical Methods](#numerical-methods)

  * [Euler's Method](#eulers-method)
  * [Fourth-Order Runge-Kutta](#fourth-order-runge-kutta)
  * [Velocity Verlet](#velocity-verlet)
  * [Convergence and Energy Conservation](#convergence-and-energy-conservation)
* [Results](#results)

  * [1. Numerical Convergence](#1-numerical-convergence)
  * [2. Long-Term Energy Conservation](#2-long-term-energy-conservation)
  * [3. Relativistic Orbital Dynamics](#3-relativistic-orbital-dynamics)
  * [4. Relativistic Precession Coefficients](#4-relativistic-precession-coefficients)
  * [5. Critical Perturbation Near the ISCO](#5-critical-perturbation-near-the-isco)
* [Key Findings](#key-findings)
* [Future Improvements](#future-improvements)
* [Conclusion](#conclusion)
* [Project Structure](#project-structure)
* [Requirements](#requirements)
* [Running the Project](#running-the-project)

## Motivation

Orbital dynamics is a useful system for studying numerical methods and the behavior of dynamical systems. For a Newtonian gravitational orbit, there is a known analytical solution in terms of Keplerian motion, which is useful for comparing numerical solutions with an exact result.

I first compare Euler's Method, RK4, and Velocity Verlet using two different orbital shapes: an eccentric orbit and a nearly circular orbit. The analytical Keplerian solution provides a reference for the numerical error. I measure this error as the timestep is changed. I then compare the long-term energy behavior of the different methods over many orbital periods.

Next, I introduce a relativistic correction to study periapsis precession. I measure the precession per orbit and fit the results to determine the leading- and next-to-leading-order coefficients. These coefficients can then be used to predict the precession for a given orbit.

Finally, I investigate the stability of circular orbits near the ISCO. Small radial perturbations are applied to the circular orbits. The perturbation size is increased until the orbit becomes unstable. This allows the numerical ISCO to be compared with the analytical prediction.

## Mathematical Model

The simulations use normalized units with

$$
GM=1.
$$

The state vector is

$$
\mathbf s=(x,y,v_x,v_y),
$$

with

$$
r=\sqrt{x^2+y^2}.
$$

### Newtonian Orbital Dynamics

Newtonian gravity gives

$$
\mathbf a=-\frac{\mathbf r}{r^3},
$$

or

$$
\ddot{x}=-\frac{x}{r^3},
\qquad
\ddot{y}=-\frac{y}{r^3}.
$$

The specific energy and angular momentum are

$$
E=\frac{v^2}{2}-\frac1r,
\qquad
L_z=xv_y-yv_x.
$$

For bound orbits, the semimajor axis and eccentricity are

$$
a=-\frac{1}{2E},
\qquad
e=\sqrt{1+2EL_z^2}.
$$

These quantities are used to construct the analytical Keplerian reference solution.

### Keplerian Analytical Solution

The Newtonian orbit is obtained using the eccentric anomaly $E_{\rm anom}$, which satisfies Kepler's equation

$$
E_{\rm anom}-e\sin E_{\rm anom}=nt+C,
$$

where $C$ encodes the initial conditions and

$$
n=\sqrt{\frac{1}{a^3}}.
$$

The position is then

$$
x=a(e-\cos E_{\rm anom}),
$$

$$
y=-a\sqrt{1-e^2}\sin E_{\rm anom}.
$$

Kepler's equation is solved numerically with `fsolve`, providing an independent reference solution for testing the numerical integrators.

### Relativistic Correction

The relativistic model adds the Schwarzschild correction to the Newtonian force:

$$
\mathbf a=
-\frac{\mathbf r}{r^3}
-\frac{3L^2\mathbf r}{r^5}.
$$

Introducing $u=1/r$, the corresponding Binet equation is

$$
\boxed{
u''+u=\frac1{L^2}+3u^2
}.
$$

This equation provides the analytical basis for the precession and ISCO analyses.

### Circular Orbits, Stability, and the ISCO

For a circular orbit, $u=u_0$ and $u''=0$, giving

$$
3u_0^2-u_0+\frac1{L^2}=0.
$$

The stable circular-orbit branch is

$$
u_0=
\frac{1-\sqrt{1-12/L^2}}{6}.
$$

To determine its stability, introduce a small radial perturbation,

$$
u=u_0+\delta u.
$$

Keeping only first-order terms gives

$$
\delta u''+(1-6u_0)\delta u=0,
$$

so the radial frequency is

$$
\omega_r^2=1-6u_0
=1-\frac6{r_c}.
$$

For $\omega_r^2>0$,

$$
\delta u\propto\cos(\omega_r\phi),
$$

so the perturbation remains bounded and the orbit is stable. For $\omega_r^2<0$, the frequency is imaginary and

$$
\delta u\propto e^{\pm\gamma\phi},
$$

so the perturbation grows and the orbit is unstable. At $\omega_r=0$, the radial restoring force vanishes, giving marginal stability.

Thus,

$$
r_c>6\Rightarrow\text{stable},
\qquad
r_c=6\Rightarrow\text{marginally stable},
\qquad
r_c<6\Rightarrow\text{unstable}.
$$

Therefore,

$$
\boxed{r_{\rm ISCO}=6}.
$$

To initialize the circular orbits used in the simulations, the circularity condition is written in terms of the tangential velocity. The required centripetal acceleration satisfies

$$
\frac{V_t^2}{r}
=\frac1{r^2}+\frac{3L^2}{r^4}.
$$

Substituting $L=rV_t$ gives

$$
V_t^2r=1+3V_t^2,
$$

so

$$
\boxed{
r_c=\frac1{V_t^2}+3
}.
$$

The corresponding initial state is

$$
(x_0,y_0,v_{x0},v_{y0})
=(r_c,0,0,V_t).
$$

The simulations then test stability by applying small radial perturbations,

$$
r_0\rightarrow r_0(1\pm\epsilon),
$$

and measuring the resulting radial deviation.

### Relativistic Precession

Since the radial perturbation oscillates as

$$
\delta u\propto\cos(\omega_r\phi),
$$

one complete radial cycle requires

$$
\Phi=\frac{2\pi}{\omega_r}
=\frac{2\pi}{\sqrt{1-6u_0}}.
$$

For a Newtonian orbit, the radial cycle closes after one full revolution, so $\Phi=2\pi$. The relativistic periapsis advance is therefore

$$
\Delta\phi=
\Phi-2\pi=
2\pi
\left[
\frac{1}{\sqrt{1-6u_0}}-1
\right].
$$

From the circular-orbit equation,

$$
u_0=
\frac{1-\sqrt{1-12/L^2}}{6},
$$

so

$$
1-6u_0=
\sqrt{1-\frac{12}{L^2}}.
$$

Therefore,

$$
\Delta\phi=
2\pi
\left[
\left(1-\frac{12}{L^2}\right)^{-1/4}
-1
\right].
$$

Using the Taylor expansion

$$
(1-x)^{-1/4}=
1+\frac{x}{4}
+\frac{5x^2}{32}
+\cdots,
$$

with $x=12/L^2$, gives

$$
\Delta\phi=
2\pi
\left[
\frac{3}{L^2}
+
\frac{45}{2L^4}
+\cdots
\right].
$$

Thus,

$$
\boxed{
\Delta\phi=
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4}
+\mathcal O(L^{-6})
}.
$$

Numerically, periapsides are identified when

$$
v_r=\frac{\mathbf r\cdot\mathbf v}{r}
$$

crosses from negative to positive. The angular separation between successive periapsides gives the measured precession, which is fitted against $1/L^2$.

## Numerical Methods

Three numerical integrators are used: Euler, fourth-order Runge-Kutta (RK4), and Velocity Verlet.

### Euler's Method

Euler's Method was derived in the previous project, where its implementation and convergence were established. It advances the state using the derivative evaluated at the current time:

$$
\mathbf s_{n+1}
=\mathbf s_n+
\mathbf f(\mathbf s_n,t_n)h.
$$

It has global error

$$
O(h).
$$

Here, Euler's Method serves primarily as a baseline for comparison.

### Fourth-Order Runge-Kutta

RK4 improves on Euler's Method by evaluating the derivative at four points within each timestep. For

$$
\frac{d\mathbf s}{dt}
=\mathbf f(\mathbf s,t),
$$

the four derivative estimates are

$$
k_1
=\mathbf f(\mathbf s_n,t_n),
$$

$$
k_2
=\mathbf f
\left(
\mathbf s_n+\frac{h}{2}k_1,
t_n+\frac{h}{2}
\right),
$$

$$
k_3
=\mathbf f
\left(
\mathbf s_n+\frac{h}{2}k_2,
t_n+\frac{h}{2}
\right),
$$

$$
k_4
=\mathbf f
\left(
\mathbf s_n+hk_3,
t_n+h
\right).
$$

These are combined to produce the next state:

$$
\mathbf s_{n+1}
=\mathbf s_n+
\frac{h}{6}
(k_1+2k_2+2k_3+k_4).
$$

RK4 has global error

$$
O(h^4).
$$

### Velocity Verlet

Velocity Verlet is designed for systems where acceleration depends on position. It first advances the position using the current velocity and acceleration:

$$
\mathbf r_{n+1}
=\mathbf r_n+
\mathbf v_nh+
\frac12\mathbf a_nh^2.
$$

The acceleration is then recomputed from the new position:

$$
\mathbf a_{n+1}
=\mathbf a(\mathbf r_{n+1}),
$$

and the velocity is updated using the average of the old and new accelerations:

$$
\mathbf v_{n+1}
=\mathbf v_n+
\frac12
(\mathbf a_n+\mathbf a_{n+1})h.
$$

Velocity Verlet has global error

$$
O(h^2),
$$

and its time-symmetric structure gives good long-term energy behavior for conservative orbital systems.

### Convergence and Energy Conservation

Numerical convergence is tested by comparing the numerical position with the Keplerian reference using the relative position error:

$$
\epsilon_r=
\frac{
\left|
\mathbf r_{\rm numerical}-
\mathbf r_{\rm exact}
\right|
}{
\left|
\mathbf r_{\rm exact}
\right|
}.
$$

The expected scaling is

$$
\epsilon_r\propto h^p,
$$

where $p$ is the method's order.

Long-term stability is measured using the relative energy error,

$$
\frac{E(t)-E(0)}{|E(0)|}\times100\%.
$$

This allows the project to compare both short-term accuracy and long-term energy behavior of the three integrators.


## Results

### 1. Numerical Convergence

![Numerical convergence](images/numerical_convergence-2.png?raw=true)

**Figure 1.** *Numerical convergence for eccentric and near-circular Newtonian orbits. The top row shows the position error on a linear scale, while the bottom row shows the same error on a log-log scale with fitted convergence slopes.*

The eccentric ($e=0.730$) and near-circular ($e=0.093$) initial conditions provide two distinct tests of the numerical integrators. The near-circular orbit tests conditions similar to the near-circular relativistic orbits studied later in the project. The eccentric orbit acts as a stress test, since the radius and velocity change much more rapidly near periapsis.

The linear-scale plots show only Velocity Verlet and RK4 because Euler's Method has much larger errors over this timestep range. This difference is especially clear in the log-log plots, where Euler's error is on the order of $10^2$. Over the short interval of a few orbital periods, RK4 clearly outperforms Velocity Verlet and maintains very small position errors, particularly for timesteps between $10^{-3}$ and $10^{-2}$.

The log-log plots also show the expected convergence behavior,

$$
\epsilon_r\propto h^p,
$$

where $p$ is the convergence order. Euler's Method shows first-order behavior, while Velocity Verlet follows the expected second-order scaling. RK4 shows much steeper convergence and reaches extremely small errors at smaller timesteps.

Overall, this experiment shows that RK4 is the most accurate method over the short timescale considered, while the eccentric orbit provides a more demanding test of the integrators than the near-circular case.

### 2. Long-Term Energy Conservation

![Energy conservation](images/energy_conservation.png?raw=true)

**Figure 2.** *Relative energy error for Euler's Method, RK4, and Velocity Verlet during orbital evolution. The simulations compare an eccentric orbit and a near-circular orbit.*

This experiment continues the comparison from the previous section by examining how the three integrators behave over many orbital periods. The eccentric ($e=0.730$) and near-circular ($e=0.093$) orbits provide two different tests of long-term energy conservation.

The simulations use $h=0.1$, with Euler and RK4 run for $500$ orbital periods and Velocity Verlet run for $10$ periods. Euler's Method shows significant energy drift in both cases, with the error becoming especially large for the eccentric orbit. The rapid changes near periapsis make the eccentric orbit more demanding and amplify the energy error over time.

RK4 performs much better, particularly for the near-circular orbit, where the energy error remains very small for most of the simulation. However, because RK4 is not symplectic, its energy error gradually accumulates over long times. This becomes much more noticeable for the eccentric orbit over $500$ periods.

Velocity Verlet shows a different behavior. Instead of steadily drifting, its energy error remains bounded and oscillates around the correct value. This gives Velocity Verlet better long-term energy conservation, especially for the eccentric orbit.

Overall, RK4 outperforms Velocity Verlet for the shorter simulations used in the convergence test, except when highly eccentric orbits are evolved for a long time. RK4 is therefore chosen as the main integrator for the relativistic simulations, which focus on near-circular orbits over relatively short timescales. Its high accuracy is important for measuring small effects such as periapsis precession.

### 3. Relativistic Orbital Dynamics

![Relativistic orbital dynamics](images/relativistic_orbits.png?raw=true)

**Figure 3.** *Relativistic orbital dynamics for weak-field, intermediate-radius, and near-ISCO orbits. Each panel compares outward and inward radial perturbations.*

The first panel shows the weak-field limit, where the relativistic correction is small and the orbit remains close to Newtonian behavior.

The second panel uses a smaller orbital radius, making the relativistic correction more important. The orbit no longer closes after a single revolution, producing a visible precession of the periapsis.

The final panel approaches the ISCO region. Here, the stability of the circular orbit becomes much more sensitive to radial perturbations. Outward and inward perturbations can produce qualitatively different behavior as the orbit approaches the stability boundary.

The simulations therefore show how the relativistic correction becomes increasingly important as the orbital radius decreases.

### 4. Relativistic Precession Coefficients

![Relativistic precession](images/relativistic_precession.png?raw=true)

**Figure 4.** *Numerical extraction of the leading-order and next-to-leading-order relativistic precession coefficients, together with a comparison between numerical precession and the theoretical expansion.*

The periapsis advance is measured by finding successive periapsides and calculating their angular separation.

The measured precession is plotted against

$$
\frac{1}{L^2}.
$$

Multiplying the measured precession by $L^2$ allows the leading-order coefficient to be extracted using a linear fit.

The first fitted coefficient approaches the theoretical value

$$
6\pi.
$$

The next-order behavior can then be examined by fitting the remaining dependence and comparing it with the predicted coefficient

$$
45\pi.
$$

The final panel compares the numerical results for $\epsilon=0.01$ with the leading-order approximation

$$
\Delta\phi=
\frac{6\pi}{L^2}
$$

and the next-to-leading-order approximation

$$
\Delta\phi=
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4}.
$$

The comparison shows how adding the higher-order correction improves the theoretical description of the numerical results away from the strict weak-field limit.

### 5. Critical Perturbation Near the ISCO

![Critical perturbation near the ISCO](images/isco_critical_perturbation.png?raw=true)

**Figure 5.** *Critical fractional radial perturbation required to produce a greater than 10% radial deviation for circular orbits near the ISCO.*

The final simulation examines the stability of circular orbits as their radius approaches the theoretical ISCO at

$$
r_c=6.
$$

For each circular-orbit radius, both inward and outward radial perturbations are tested. The perturbation size is increased until the orbit's radial deviation exceeds 10%.

The simulation is repeated for different maximum integration times, measured in orbital periods.

The resulting curves show how the critical perturbation changes as the circular orbit approaches the ISCO. The decreasing stability margin near the theoretical boundary provides a numerical indication of the transition from stable to unstable circular motion.

Using different integration times is also important because an unstable orbit may require time for the perturbation to grow enough to cross the chosen 10% threshold.

## Key Findings

* Euler's Method is first-order accurate and shows significant long-term energy drift in orbital simulations.
* RK4 achieves fourth-order convergence and maintains much smaller numerical errors over long integrations.
* Velocity Verlet is second-order accurate but has favorable long-term behavior for conservative orbital systems.
* The Newtonian simulations can be compared directly with an analytical Keplerian solution obtained from Kepler's equation.
* The relativistic correction causes the periapsis to advance, so the orbit is no longer exactly closed.
* The leading-order numerical precession coefficient approaches the theoretical value $6\pi$.
* The next-to-leading-order correction is consistent with a $45\pi/L^4$ contribution.
* Relativistic effects become increasingly important as the orbital radius decreases.
* Circular orbits become increasingly sensitive to radial perturbations as the radius approaches the ISCO at $r_c=6$.

## Future Improvements

* Compare the numerical trajectories with a higher-precision reference solution rather than relying only on the analytical Newtonian solution.
* Implement adaptive timestep integration for the relativistic simulations.
* Compare RK4 and Velocity Verlet over the same number of orbital periods for a more direct long-term comparison.
* Investigate angular-momentum conservation in addition to energy conservation.
* Improve the periapsis detection by interpolating the radial-velocity zero crossing more accurately.
* Explore the dependence of relativistic precession on eccentricity as well as angular momentum.
* Extend the relativistic model to a more complete treatment of Schwarzschild geodesics.
* Investigate the growth rate of perturbations near the ISCO rather than only measuring whether a fixed 10% threshold is crossed.

## Conclusion

This project uses orbital dynamics to study numerical integration, conservation laws, relativistic corrections, perihelion precession, and orbital stability.

The Newtonian problem provides a useful benchmark because its analytical Keplerian solution allows numerical errors to be measured directly. The convergence analysis shows the different orders of accuracy of Euler's Method, RK4, and Velocity Verlet, while the long-term simulations demonstrate that accuracy and physical stability are not necessarily the same thing.

The relativistic simulations then show how a small modification to the effective potential changes the qualitative behavior of the orbit. Instead of closing exactly after each revolution, the periapsis advances. Measuring this angular shift allows the leading-order and next-to-leading-order coefficients of the precession to be extracted numerically.

The final simulations investigate circular-orbit stability near the ISCO. As the circular radius approaches the theoretical value $r_c=6$, small radial perturbations become increasingly important, providing a numerical picture of the transition toward unstable orbital motion.

Overall, the project connects numerical methods with physical predictions. The same simulation framework is used first to test the accuracy and long-term behavior of numerical integrators and then to investigate phenomena that have no simple closed Newtonian orbit, such as relativistic precession and the stability boundary near the ISCO.

## Project Structure

```text
Orbital-Dynamics-and-Relativistic-Precession/
├── main.py
├── physics.py
├── solvers.py
├── analysis.py
├── numerical_convergence.png
├── energy_conservation.png
├── relativistic_orbits.png
├── relativistic_precession.png
├── isco_critical_perturbation.png
├── requirements.txt
└── README.md
```

## Requirements

* Python 3
* NumPy
* Matplotlib
* SciPy

## Running the Project

### 1. Clone the repository

```bash
git clone <repository-url>
cd <repository-directory>
```

### 2. Install the dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the project

```bash
python main.py
```
