# Orbital Dynamics and Relativistic Precession

A computational physics project using numerical integration to study Newtonian orbital motion, numerical convergence, long-term energy conservation, relativistic orbital precession, and stability near the innermost stable circular orbit (ISCO).

The project compares Euler's Method, fourth-order Runge-Kutta (RK4), and Velocity Verlet for Newtonian orbital dynamics. Analytical Keplerian orbits are used as a reference solution to measure numerical error and convergence. The project then introduces a relativistic correction to the gravitational potential and uses numerical simulations to investigate perihelion precession and the stability of circular orbits near the ISCO.

![Numerical convergence](images/numerical_convergence.png?raw=true)

**Figure 1.** *Numerical convergence of the orbital integrators for eccentric and near-circular Newtonian orbits.*

## Table of Contents

* [Motivation](#motivation)
* [Mathematical Model](#mathematical-model)

  * [Newtonian Orbital Dynamics](#newtonian-orbital-dynamics)
  * [Orbital Energy and Angular Momentum](#orbital-energy-and-angular-momentum)
  * [Keplerian Analytical Solution](#keplerian-analytical-solution)
  * [Relativistic Correction](#relativistic-correction)
  * [Circular Orbits and the ISCO](#circular-orbits-and-the-isco)
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

Orbital dynamics provides a useful way to study both numerical integration and the underlying physics of dynamical systems.

For a Newtonian gravitational orbit, the equations of motion have a known analytical solution in terms of Keplerian motion. This makes the system particularly useful for testing numerical methods because the numerical solution can be compared directly with an exact result.

I first compare Euler's Method, RK4, and Velocity Verlet using two different orbital shapes: an eccentric orbit and a nearly circular orbit. The analytical Keplerian solution provides a reference position, allowing the numerical error to be measured as the timestep is changed.

I then examine the long-term energy behavior of the different integrators. This is especially important for orbital simulations because a method can have small short-term errors while still producing unphysical behavior over many orbital periods.

Finally, I introduce a relativistic correction to the effective potential. This causes the orbit to precess rather than closing on itself. By measuring the angular shift between successive periapsides, the numerical simulations can be used to extract the leading-order and next-to-leading-order relativistic precession coefficients.

The final part of the project investigates the stability of circular orbits near the ISCO. Small radial perturbations are applied to circular orbits to determine when the orbit becomes unstable.

## Mathematical Model

The simulations use normalized units in which the gravitational parameter is

$$
GM=1.
$$

The position and velocity of the orbiting body are written as

$$
\mathbf r=
\begin{bmatrix}
x\\
y
\end{bmatrix},
\qquad
\mathbf v=
\begin{bmatrix}
v_x\\
v_y
\end{bmatrix}.
$$

The full numerical state is therefore

$$
\mathbf s=
\begin{bmatrix}
x\\
y\\
v_x\\
v_y
\end{bmatrix}.
$$

The distance from the central object is

$$
r=\sqrt{x^2+y^2}.
$$

### Newtonian Orbital Dynamics

For Newtonian gravity, the acceleration is directed toward the origin and has magnitude proportional to $1/r^2$:

$$
\mathbf a=-\frac{GM}{r^3}\mathbf r.
$$

With $GM=1$, this becomes

$$
\mathbf a=-\frac{\mathbf r}{r^3}.
$$

In Cartesian coordinates,

$$
\ddot{x}=-\frac{x}{r^3},
$$

$$
\ddot{y}=-\frac{y}{r^3}.
$$

The corresponding first-order system is

$$
\frac{dx}{dt}=v_x,
$$

$$
\frac{dy}{dt}=v_y,
$$

$$
\frac{dv_x}{dt}=-\frac{x}{r^3},
$$

$$
\frac{dv_y}{dt}=-\frac{y}{r^3}.
$$

These equations are implemented by `orbit_derivative()` in `physics.py`.

For the Newtonian simulations, the initial position is

$$
\mathbf r_0=(3,0),
$$

while the initial velocity is chosen as

$$
\mathbf v_0=(0,V_{y0}).
$$

Two values are used:

$$
V_{y0}=0.3
$$

for the eccentric orbit and

$$
V_{y0}=0.55
$$

for the near-circular orbit.

### Orbital Energy and Angular Momentum

For a unit-mass body, the specific mechanical energy is

$$
E=\frac{v^2}{2}-\frac{GM}{r}.
$$

With the normalized value $GM=1$,

$$
E=\frac{v^2}{2}-\frac1r.
$$

For a bound orbit, $E<0$.

The specific angular momentum is

$$
L_z=xv_y-yv_x.
$$

Because the Newtonian gravitational force is central, angular momentum is conserved.

The eccentricity can be calculated from the energy and angular momentum:

$$
e=
\sqrt{
1+\frac{2EL_z^2}{(GM)^2}
}.
$$

This quantity determines the shape of the Keplerian orbit.

For the initial conditions used in the simulations, the $V_{y0}=0.3$ orbit is significantly eccentric, while the $V_{y0}=0.55$ orbit is close to circular.

For a bound Newtonian orbit, the semi-major axis is related to the energy by

$$
a=-\frac{GM}{2E}.
$$

The orbital period is then

$$
T=2\pi\sqrt{\frac{a^3}{GM}}.
$$

This is the form of Kepler's third law used to determine the simulation times for the long-term energy comparison.

### Keplerian Analytical Solution

The Newtonian two-body problem has an analytical solution in terms of the eccentric anomaly $E_{\mathrm{anom}}$.

The mean motion is

$$
n=\sqrt{\frac{GM}{a^3}},
$$

and Kepler's equation is

$$
E_{\mathrm{anom}}-e\sin E_{\mathrm{anom}}=
nt+C,
$$

where $C$ is determined by the initial conditions.

Once the eccentric anomaly has been found, the position can be written as

$$
x=a(e-\cos E_{\mathrm{anom}})
$$

and

$$
y=-a\sqrt{1-e^2}\sin E_{\mathrm{anom}}.
$$

The project solves Kepler's equation numerically using `fsolve` and uses the resulting position as the analytical reference.

This gives an independent solution against which the numerical integrators can be tested.

### Relativistic Correction

To investigate relativistic effects, the Newtonian potential is modified by an additional angular-momentum-dependent term.

The effective potential used in the project is

$$
V(r)=-\frac1r+\frac{L^2}{r^3},
$$

where $L$ is the angular momentum associated with the orbit.

The radial derivative of this potential gives the corresponding acceleration:

$$
\mathbf a=
-\frac{\mathbf r}{r^3}-
\frac{3L^2\mathbf r}{r^5}.
$$

In Cartesian coordinates,

$$
\ddot{x}=
-\frac{x}{r^3}-
\frac{3L^2x}{r^5},
$$

and

$$
\ddot{y}=
-\frac{y}{r^3}-
\frac{3L^2y}{r^5}.
$$

The resulting first-order equations are therefore

$$
\frac{dx}{dt}=v_x,
$$

$$
\frac{dy}{dt}=v_y,
$$

$$
\frac{dv_x}{dt}=
-\frac{x}{r^3}-
\frac{3L^2x}{r^5},
$$

$$
\frac{dv_y}{dt}=
-\frac{y}{r^3}-
\frac{3L^2y}{r^5}.
$$

This model is implemented by `relativistic_orbit_derivative()`.

Unlike the Newtonian problem, the additional term changes the orbital dynamics so that the periapsis does not return to exactly the same angular position after each orbit. The resulting advance of the periapsis is the relativistic precession measured in the later simulations.

### Circular Orbits and the ISCO

Circular orbits are particularly useful for studying orbital stability.

For the relativistic model, the initial circular-orbit radius is chosen as

$$
r_c=\frac{1}{V_{y0}^2}+3.
$$

The corresponding initial state is

$$
(x_0,y_0,v_{x0},v_{y0})=
(r_c,0,0,V_{y0}).
$$

As the orbital velocity increases, the circular-orbit radius decreases. This allows the simulations to approach the region where circular orbits become unstable.

The theoretical ISCO in the normalized units used here occurs at

$$
r_c=6.
$$

To investigate this transition numerically, the circular orbit is perturbed radially by a small fractional amount:

$$
r_0\rightarrow r_0(1\pm\epsilon).
$$

Both outward and inward perturbations are tested.

If the orbit remains close to the circular trajectory, the perturbation is stable. If the radial deviation grows substantially, the circular orbit is unstable.

The final simulation measures the smallest perturbation that produces a radial deviation greater than 10% over a specified number of orbital periods.

### Relativistic Precession

A bound Newtonian orbit is closed after one orbital period. The relativistic correction breaks this exact closure, causing the periapsis to advance from one orbit to the next.

The periapsis angle is found by detecting when the radial velocity changes from negative to positive:

$$
v_r=\frac{\mathbf r\cdot\mathbf v}{r}.
$$

At each periapsis, the angular position is calculated using

$$
\phi=\tan^{-1}\left(\frac{y}{x}\right).
$$

The angles are unwrapped so that the total accumulated angular motion can be measured continuously.

The average angular separation between successive periapsides is then calculated. This gives the numerical periapsis-to-periapsis angular advance.

For the weak-field expansion used in the analysis, the precession is expected to have the form

$$
\Delta\phi
\approx
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4}
+\cdots.
$$

The simulations test these coefficients numerically by fitting the measured precession as a function of

$$
\frac{1}{L^2}.
$$

## Numerical Methods

The equations are integrated numerically using three different methods: Euler's Method, fourth-order Runge-Kutta, and Velocity Verlet.

Comparing these methods allows both their formal accuracy and their long-term physical behavior to be studied.

### Euler's Method

Euler's Method is the simplest of the numerical integrators used in the project.

For a general state vector,

$$
\mathbf s_{n+1}=
\mathbf s_n+
\mathbf f(\mathbf s_n,t_n)\Delta t.
$$

For the orbital equations, this gives

$$
x_{n+1}=x_n+v_{x,n}\Delta t,
$$

$$
y_{n+1}=y_n+v_{y,n}\Delta t,
$$

$$
v_{x,n+1}=
v_{x,n}
-\frac{x_n}{r_n^3}\Delta t,
$$

and

$$
v_{y,n+1}=
v_{y,n}
-\frac{y_n}{r_n^3}\Delta t.
$$

Euler's Method has local truncation error

$$
O(\Delta t^2)
$$

and global error

$$
O(\Delta t).
$$

It is computationally inexpensive, but its errors accumulate over time. For orbital dynamics, this can lead to significant changes in the orbital energy and eventually distort the orbit.

### Fourth-Order Runge-Kutta

The fourth-order Runge-Kutta method, or RK4, uses four derivative evaluations per timestep.

The four intermediate slopes are

$$
k_1=f(s_n,t_n),
$$

$$
k_2=f\left(s_n+\frac{\Delta t}{2}k_1,t_n+\frac{\Delta t}{2}\right),
$$

$$
k_3=f\left(s_n+\frac{\Delta t}{2}k_2,t_n+\frac{\Delta t}{2}\right),
$$

and

$$
k_4=f(s_n+\Delta t k_3,t_n+\Delta t).
$$

The final update is

$$
s_{n+1}=
s_n+
\frac{\Delta t}{6}
(k_1+2k_2+2k_3+k_4).
$$

RK4 has local truncation error

$$
O(\Delta t^5)
$$

and global error

$$
O(\Delta t^4).
$$

This makes it substantially more accurate than Euler's Method for sufficiently small timesteps.

RK4 is used throughout the relativistic simulations because of its high accuracy and general applicability to the modified equations of motion.

### Velocity Verlet

Velocity Verlet is a second-order integration method that updates the position using the current acceleration and then uses the new acceleration to update the velocity.

The position update is

$$
\mathbf r_{n+1}=
\mathbf r_n+
\mathbf v_n\Delta t+
\frac12\mathbf a_n\Delta t^2.
$$

The acceleration is then recalculated at the new position,

$$
\mathbf a_{n+1}=
\mathbf a(\mathbf r_{n+1}),
$$

and the velocity is updated using

$$
\mathbf v_{n+1}=
\mathbf v_n+
\frac12
(\mathbf a_n+\mathbf a_{n+1})
\Delta t.
$$

Velocity Verlet has global error

$$
O(\Delta t^2).
$$

Its main advantage for this project is its behavior for conservative orbital systems. Although it has lower formal order than RK4, it has favorable long-term energy behavior and does not require the acceleration to be written as a first-order state derivative.

This makes it particularly useful for comparing short-term numerical accuracy with long-term physical stability.

### Convergence and Energy Conservation

Two different properties of the numerical methods are tested.

The first is **convergence**. The final numerical position is compared with the analytical Keplerian position for a range of timestep sizes,

$$
0.002\leq\Delta t\leq0.1.
$$

The position error is

$$
\epsilon_r=
\left|
\mathbf r_{\mathrm{numerical}}-
\mathbf r_{\mathrm{exact}}
\right|.
$$

A log-log plot of error against timestep should approximately follow

$$
\epsilon_r\propto\Delta t^p,
$$

where $p$ is the order of convergence.

The second property is **long-term energy conservation**.

The relative energy error is calculated as

$$
\frac{E(t)-E(0)}{|E(0)|}\times100\%.
$$

The simulations compare Euler and RK4 over 500 orbital periods and Velocity Verlet over a shorter interval.

This distinction is useful because a method can have excellent local accuracy while still displaying undesirable long-term behavior.

## Results

### 1. Numerical Convergence

![Numerical convergence](images/numerical_convergence.png?raw=true)

**Figure 1.** *Numerical convergence for eccentric and near-circular Newtonian orbits. The top row shows the position error on a linear scale, while the bottom row shows the same error on a log-log scale with fitted convergence slopes.*

The eccentric and near-circular initial conditions provide two different tests of the numerical methods.

The log-log plots make the convergence rate particularly clear. If the numerical error follows

$$
\epsilon_r\propto\Delta t^p,
$$

then the slope of the log-log curve approaches $p$.

Euler's Method is expected to show first-order convergence, while RK4 should show fourth-order convergence. Velocity Verlet is expected to show second-order convergence.

The comparison also demonstrates that convergence is not only a property of the numerical method but can depend on the dynamics being simulated. The eccentric orbit contains stronger variations in velocity and radius, making it a more demanding numerical problem than the nearly circular orbit.

### 2. Long-Term Energy Conservation

![Energy conservation](images/energy_conservation.png?raw=true)

**Figure 2.** *Relative energy error for Euler's Method, RK4, and Velocity Verlet during orbital evolution. The simulations compare an eccentric orbit and a near-circular orbit.*

The energy plots show how numerical errors accumulate over many orbital periods.

Euler's Method generally exhibits significant long-term energy drift. Because the orbit is repeatedly updated using only the derivative at the beginning of each timestep, small errors accumulate and can alter the orbital energy substantially.

RK4 has much smaller local errors and therefore maintains the orbital energy much more accurately over the same interval.

Velocity Verlet provides a different type of advantage. Its second-order accuracy is lower than RK4, but its structure is well suited to conservative systems. Its energy error remains bounded rather than growing in the same way as a simple first-order method.

The difference between the eccentric and near-circular cases is also useful. Eccentric orbits experience larger variations in radius and velocity, so numerical errors can have a more pronounced effect on the dynamics.

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
