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

Orbital dynamics is a useful system for studying numerical methods and the behavior of dynamical systems. For a Newtonian gravitational orbit, there is a known analytical solution in terms of Keplerian motion, which is useful for comparing numerical solutions to an exact result.

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

The Binet equation

$$
u''+u=\frac1{L^2}+3u^2
$$

is nonlinear because of the $3u^2$ term. In the weak-field limit, $u=1/r$ is small, so the relativistic correction is also small. This suggests treating the relativistic effects as a perturbation of the Newtonian orbit and expanding the resulting precession in powers of the small parameter $1/L^2$. The precession can therefore be written as

$$
\frac{c_1}{L^2}
+
\frac{c_2}{L^4}
+
\mathcal O(L^{-6}),
$$

where $c_1$ is the leading-order relativistic correction and $c_2$ is the next-to-leading-order correction.

For a circular orbit, the radial frequency is

$$
\omega_r^2=1-6u_0,
$$

so one complete radial cycle requires

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

Since the weak-field limit corresponds to $1/L^2\ll1$, this expression can be expanded as a power series in $1/L^2$:

$$
(1-x)^{-1/4}=
1+\frac{x}{4}
+\frac{5x^2}{32}
+\cdots,
$$

with $x=12/L^2$. This gives

$$
\Delta\phi=
2\pi
\left[
\frac{3}{L^2}
+
\frac{45}{2L^4}
+\cdots
\right],
$$

and therefore

$$
\boxed{
\Delta\phi=
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4}
+\mathcal O(L^{-6})
}.
$$

This expansion gives a natural way to test the numerical results: the coefficients of the $1/L^2$ and $1/L^4$ terms can be extracted from the simulations and compared with the theoretical predictions.

Numerically, periapsides are identified when

$$
v_r=\frac{\mathbf r\cdot\mathbf v}{r}
$$

crosses from negative to positive. The angular separation between successive periapsides gives the measured precession, which is then fitted against $1/L^2$.

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

The expected scaling is $\epsilon_r\propto h^p,$ where $p$ is the method's order.

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

This experiment examines relativistic orbital behavior at three different radii, using both outward ($+\epsilon$) and inward ($-\epsilon$) radial perturbations.

In the weak-field case, with $V_{y0}=0.01$ and an orbital radius of approximately $r=104$, relativistic corrections are very small. Both perturbations remain stable and produce closed, Newtonian-like elliptical orbits with little visible precession.

At the intermediate radius, using $V_{y0}=0.15$ and $r\approx50$, relativistic effects become more noticeable. Both perturbations remain bounded, but the orbits no longer close after one revolution. Instead, the periapsis gradually shifts, producing the characteristic rosette pattern associated with relativistic precession.

Near the ISCO, with $V_{y0}=0.6$ and $r\approx6$, the behavior changes significantly. The outward perturbation moves the orbit into the stable region and produces a bounded, precessing trajectory. In contrast, the inward perturbation pushes the particle into the unstable region, causing it to plunge inward toward the central mass.

The ISCO plot also shows a numerical artifact in the inward trajectory. The trajectory is terminated near the origin because the relativistic force becomes extremely large at small $r$. With a fixed timestep, the integrator cannot accurately resolve this steep force gradient, causing the numerical solution to behave incorrectly near $r=0$ rather than smoothly approaching the central singularity.

### 4. Relativistic Precession Coefficients

![Relativistic precession coefficients](images/precession_coefficients.png?raw=true)

**Figure 4.** *Numerical extraction of the leading-order and next-to-leading-order relativistic precession coefficients.*

Figure 4 shows how the coefficients in the relativistic precession expansion are extracted numerically. Starting from

$$
\Delta\phi=
\frac{c_1}{L^2}
+
\frac{c_2}{L^4}
+
\mathcal O(L^{-6}),
$$

the equation is multiplied by $L^2$ and written as

$$
y=c_1+c_2x+O(x^{2}),
$$

where

$$
y=\Delta\phi L^2,
\qquad
x=\frac1{L^2}.
$$

A linear fit is then performed using between 3 and 9 of the weakest-field points closest to $x=0$, corresponding to $L\to\infty$. The intercept gives the leading-order coefficient $c_1$, while the slope gives the next-to-leading-order coefficient $c_2$.

The extracted $c_1$ values converge very closely to the theoretical prediction

$$
c_1=6\pi\approx18.85,
$$

and remain consistent across the tested eccentricities, $\epsilon=0.01$ to $0.10$. The extraction of $c_2$ is more sensitive to higher-order terms because points farther from the weak-field limit contain larger $\mathcal O(L^{-6})$ contributions. Using fewer points closer to $x=0$ gives values that cluster around the theoretical prediction

$$
c_2=45\pi\approx141.37.
$$

![Numerical vs theoretical precession](images/numerical_vs_theoretical_precession-2.png?raw=true)

**Figure 5.** *Comparison between the numerical periapsis precession and the theoretical relativistic expansion.*

Figure 5 compares the numerical periapsis precession with the theoretical expansion. The leading-order approximation,

$$
\Delta\phi=\frac{6\pi}{L^2},
$$

agrees closely with the numerical results in the weak-field limit, where $1/L^2$ is small. As $1/L^2$ increases, however, it increasingly underestimates the measured precession.

Including the next-to-leading-order term gives

$$
\Delta\phi=
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4}.
$$

This produces much better agreement with the numerical data across the plotted range and captures the upward curvature that becomes more apparent as the gravitational field strengthens.

### 5. Critical Perturbation Near the ISCO

![Critical perturbation near the ISCO](images/isco_critical_perturbation-2.png?raw=true)

**Figure 5.** *Critical fractional radial perturbation required to produce a greater than 10% radial deviation for circular orbits near the ISCO.*

This experiment examines the numerical stability boundary near the theoretical ISCO at $r_c=6$. A $10%$ relative radial deviation is used as the numerical criterion for stability. This threshold is arbitrary and does not represent a physical boundary. The critical perturbation $\epsilon_{\rm crit}$ is defined as the largest inward fractional perturbation that stays below this threshold.

Only inward perturbations are considered because outward perturbations move the particle into a more stable region, where the orbit remains bounded around a larger radius. Perturbations with $\epsilon<\epsilon_{\rm crit}$ remain within the $10%$ criterion, while larger perturbations can push the orbit into the unstable region and cause a plunge.

The stability boundary is tested for maximum integration times of $10T$, $20T$, $30T$, and $40T$, where $T$ is the orbital period. At shorter times, slowly growing instabilities near the ISCO may not have enough time to reach the $10%$ threshold, making the stability region appear larger.

The $20T$, $30T$, and $40T$ curves converge closely, showing that the numerical stability boundary becomes largely independent of integration time. This also shows that the specific $10%$ threshold is not important once the integration time is long enough for the instability to develop.

The converged boundary drops sharply as $r_c$ approaches the ISCO and approaches zero near $r_c\approx6.1$. The small offset from the theoretical value $r_c=6$ is expected from the finite perturbation sizes and integration timestep. Overall, the results show that the critical perturbation becomes increasingly small near the ISCO, consistent with the loss of radial stability predicted by the analytical model.

## Key Findings

* Euler's Method is first-order accurate and shows significant long-term energy drift, while RK4 has fourth-order convergence and much smaller short-term errors.
* Velocity Verlet is second-order accurate but has better bounded energy behavior over long integrations of conservative orbits.
* The Newtonian simulations agree with an analytical Keplerian solution, providing a direct reference for numerical error and convergence.
* The relativistic correction causes periapsis precession, with the leading-order term consistent with $6\pi/L^2$ and the next-order term consistent with $45\pi/L^4$.
* Relativistic effects become stronger at smaller orbital radii, and circular orbits become increasingly sensitive to perturbations as $r_c$ approaches the ISCO at $r_c=6$.

## Future Improvements

* Explore how relativistic periapsis precession depends on eccentricity and angular momentum.
* Investigate zoom-whirl behavior for near-critical relativistic orbits.
* Analyze the Schwarzschild effective potential and relate its structure to the numerically observed orbital behavior.
* Investigate relativistic gravitational scattering and the transition between scattering and capture.
* Extend the model to null geodesics and investigate gravitational light deflection.
* Study the separatrix between bound, zoom-whirl, scattering, and plunging trajectories.

## Conclusion

This project uses orbital dynamics to study numerical integration, conservation laws, relativistic corrections, periapsis precession, and orbital stability.

The Newtonian simulations provide a useful benchmark through the analytical Keplerian solution. Comparing Euler's Method, RK4, and Velocity Verlet shows the difference between numerical accuracy and long-term stability, with RK4 giving the smallest short-term errors and Velocity Verlet showing better bounded energy behavior.

The relativistic simulations then show how a small correction to the Newtonian dynamics produces periapsis precession. Measuring this shift allows the leading-order and next-to-leading-order precession coefficients to be extracted numerically and compared with the theoretical values.

Finally, the stability analysis near the ISCO shows that radial perturbations become increasingly important as $r_c$ approaches $6$. Overall, the project connects numerical methods with physical predictions and shows how the same simulation framework can be used to study both numerical behavior and relativistic orbital dynamics.

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
git clone https://github.com/rohitkb7782/Orbital-Dynamics-and-Relativistic-Precession.git
cd Orbital-Dynamics-and-Relativistic-Precession
```

### 2. Install the dependencies

```bash
pip install -r requirements.txt
```

### 3. Run the project

```bash
python main.py
```
