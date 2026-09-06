# Orbital Dynamics and Relativistic Precession

A computational physics project using numerical simulation to study orbital dynamics and relativistic effects. I compare Euler's Method, fourth-order Runge-Kutta (RK4), and Velocity Verlet for Newtonian orbits, using the analytical Keplerian solution to measure numerical error and convergence. I also compare the long-term energy behavior of the different methods.

I then add a relativistic correction to study periapsis precession and the stability of circular orbits near the innermost stable circular orbit (ISCO).

![Relativistic orbital dynamics](images/relativistic_orbits.png?raw=true)

**Figure 3.** *Relativistic orbital dynamics for weak-field, intermediate-radius, and near-ISCO orbits. Each panel compares outward and inward radial perturbations.*

## Table of Contents

* [Motivation](#motivation)
* [Mathematical Model and Numerical Methods](#mathematical-model-and-numerical-methods)
* [Results](#results)
* [Key Findings](#key-findings)
* [Future Improvements](#future-improvements)
* [Conclusion](#conclusion)
* [Project Structure](#project-structure)
* [Requirements](#requirements)
* [Running the Project](#running-the-project)

## Motivation

Orbital dynamics is a useful system for comparing numerical methods because Newtonian orbits have a known analytical solution. This gives me a way to check the numerical results directly.

I start by comparing Euler's Method, RK4, and Velocity Verlet using two different orbital shapes: an eccentric orbit and a nearly circular orbit. I use the analytical Keplerian solution as a reference and measure how the numerical error changes with the timestep. I then compare the energy of each orbit over many orbital periods to see how the methods behave over longer simulations.

Next, I add a relativistic correction to study periapsis precession. I measure the precession per orbit for different angular momenta and fit the results to extract the leading- and next-to-leading-order coefficients predicted by the theory.

Finally, I investigate circular orbits near the ISCO. I apply small radial perturbations and increase their size until the orbit becomes unstable. This gives a numerical estimate of the stability boundary that can be compared with the analytical prediction.

## Mathematical Model and Numerical Methods

The simulations use normalized units with \(GM=1\), with the state vector

$$
\mathbf{s}=(x,y,v_x,v_y).
$$

The Newtonian and relativistic equations of motion are integrated using Euler's Method, RK4, and Velocity Verlet. The underlying equations and numerical methods are introduced in:

- [Mathematical Model](docs/mathematical_model.md)
- [Numerical Methods](docs/numerical_methods.md)

## Results

### 1. Numerical Convergence

![Numerical convergence](images/numerical_convergence-2.png?raw=true)

**Figure 1.** *Numerical convergence for eccentric and near-circular Newtonian orbits. The top row shows the position error on a linear scale, while the bottom row shows the same error on a log-log scale with fitted convergence slopes.*

I use an eccentric orbit ($e=0.730$) and a near-circular orbit ($e=0.093$) to test the three integrators. The near-circular case is also similar to the orbits used later in the relativistic simulations. The eccentric orbit is more demanding because the radius and velocity change much more rapidly near periapsis.

The linear-scale plots only show Velocity Verlet and RK4 because Euler's Method has much larger errors over this timestep range. This is especially clear in the log-log plots, where Euler's error is on the order of $10^2$.

Over the short interval of a few orbital periods, RK4 gives the smallest errors, particularly for timesteps between $10^{-3}$ and $10^{-2}$. Velocity Verlet follows its expected second-order behavior, while Euler shows first-order convergence.

The log-log plots also show the expected scaling,

$$
\epsilon_r\propto h^p.
$$

Overall, RK4 is the most accurate method over the short timescale used here, while the eccentric orbit provides a more demanding test than the near-circular orbit.

### 2. Long-Term Energy Conservation

![Energy conservation](images/energy_conservation.png?raw=true)

**Figure 2.** *Relative energy error for Euler's Method, RK4, and Velocity Verlet during orbital evolution. The simulations compare an eccentric orbit and a near-circular orbit.*

I next compare the energy behavior of the three methods over many orbital periods. The same eccentric ($e=0.730$) and near-circular ($e=0.093$) orbits are used.

The simulations use $h=0.1$, with Euler and RK4 run for $500$ orbital periods and Velocity Verlet run for $10$ periods. Euler's Method shows significant energy drift in both cases, with the error becoming especially large for the eccentric orbit. The rapid changes near periapsis make the eccentric orbit more difficult to integrate accurately.

RK4 performs much better, especially for the near-circular orbit, where the energy error stays very small for most of the simulation. However, RK4 is not symplectic, so its energy error can accumulate over very long integrations. This becomes more noticeable for the eccentric orbit over $500$ periods.

Velocity Verlet behaves differently. Instead of steadily drifting, its energy error remains bounded and oscillates around the correct value. This gives it better long-term energy behavior, particularly for the eccentric orbit.

For the relativistic simulations, I use RK4 as the main integrator because those simulations focus on near-circular orbits over relatively short timescales. Its higher accuracy is useful when measuring small effects such as periapsis precession.

### 3. Relativistic Orbital Dynamics

![Relativistic orbital dynamics](images/relativistic_orbits.png?raw=true)

**Figure 3.** *Relativistic orbital dynamics for weak-field, intermediate-radius, and near-ISCO orbits. Each panel compares outward and inward radial perturbations.*

I test relativistic orbital behavior at three different radii using both outward ($+\epsilon$) and inward ($-\epsilon$) radial perturbations.

In the weak-field case, with $V_{y0}=0.01$ and an orbital radius of approximately $r=104$, the relativistic correction is very small. Both perturbations remain stable and produce nearly closed, Newtonian-like elliptical orbits with little visible precession.

At the intermediate radius, using $V_{y0}=0.15$ and $r\approx50$, the relativistic effects are easier to see. Both perturbations remain bounded, but the orbits no longer close after one revolution. Instead, the periapsis gradually shifts, producing the rosette pattern expected from relativistic precession.

Near the ISCO, with $V_{y0}=0.6$ and $r\approx6$, the behavior changes significantly. The outward perturbation moves the orbit into the stable region and produces a bounded, precessing trajectory. The inward perturbation moves the particle into the unstable region, causing it to plunge inward toward the central mass.

The inward trajectory also shows a numerical artifact near the origin. As $r$ becomes very small, the relativistic force becomes extremely large. With a fixed timestep, the integrator cannot accurately resolve this steep force, so the numerical trajectory becomes unreliable near $r=0$.

### 4. Relativistic Precession Coefficients

![Relativistic precession coefficients](images/precession_coefficients.png?raw=true)

**Figure 4.** *Numerical extraction of the leading-order and next-to-leading-order relativistic precession coefficients.*

I extract the precession coefficients by starting with

$$
\Delta\phi=
\frac{c_1}{L^2}
+
\frac{c_2}{L^4}
+
\mathcal O(L^{-6}).
$$

Multiplying by $L^2$ gives

$$
y=c_1+c_2x+O(x^{2}),
$$

where

$$
y=\Delta\phi L^2,
\qquad
x=\frac1{L^2}.
$$

I then perform a linear fit using between 3 and 9 of the weakest-field points closest to $x=0$, corresponding to $L\to\infty$. The intercept gives $c_1$, while the slope gives $c_2$.

The extracted $c_1$ values converge very closely to the theoretical prediction

$$
c_1=6\pi\approx18.85,
$$

and remain consistent across the tested eccentricities, $\epsilon=0.01$ to $0.10$.

The extraction of $c_2$ is more sensitive to higher-order terms because points farther from the weak-field limit contain larger $\mathcal O(L^{-6})$ contributions. Using fewer points closer to $x=0$ gives values that cluster around

$$
c_2=45\pi\approx141.37.
$$

![Numerical vs theoretical precession](images/numerical_vs_theoretical_precession-2.png?raw=true)

**Figure 5.** *Comparison between the numerical periapsis precession and the theoretical relativistic expansion.*

I also compare the numerical precession directly with the theoretical expansion.

The leading-order approximation,

$$
\Delta\phi=\frac{6\pi}{L^2},
$$

agrees closely with the numerical results in the weak-field limit, where $1/L^2$ is small. As $1/L^2$ increases, however, it begins to underestimate the measured precession.

Including the next-to-leading-order term,

$$
\Delta\phi=
\frac{6\pi}{L^2}
+
\frac{45\pi}{L^4},
$$

gives much better agreement across the plotted range. It also captures the upward curvature that becomes more noticeable as the gravitational field gets stronger.

### 5. Critical Perturbation Near the ISCO

![Critical perturbation near the ISCO](images/isco_critical_perturbation-2.png?raw=true)

**Figure 6.** *Critical fractional radial perturbation required to produce a greater than 10% radial deviation for circular orbits near the ISCO.*

I next look at the stability boundary near the theoretical ISCO at $r_c=6$. I use a $10%$ relative radial deviation as the numerical stability criterion. This threshold is arbitrary and is only used to define a consistent numerical measure of when the perturbation has become large.

The critical perturbation $\epsilon_{\rm crit}$ is the largest inward fractional perturbation that stays below this threshold.

I only consider inward perturbations because outward perturbations move the particle into a more stable region, where the orbit remains bounded around a larger radius. For $\epsilon<\epsilon_{\rm crit}$, the orbit stays within the $10%$ criterion. Larger perturbations can push the orbit into the unstable region and lead to a plunge.

I test maximum integration times of $10T$, $20T$, $30T$, and $40T$, where $T$ is the orbital period. At shorter integration times, slowly growing instabilities near the ISCO may not have enough time to reach the $10%$ threshold, making the stability region appear larger.

The $20T$, $30T$, and $40T$ curves converge closely, showing that the numerical stability boundary becomes mostly independent of integration time once the simulation is long enough.

The converged boundary drops sharply as $r_c$ approaches the ISCO and approaches zero near $r_c\approx6.1$. The small offset from the theoretical value $r_c=6$ is expected because of the finite perturbation sizes and timestep.

Overall, the results show that the critical perturbation becomes smaller as the orbit approaches the ISCO, consistent with the loss of radial stability predicted by the analytical model.

## Key Findings

* Euler's Method is first-order accurate and shows significant long-term energy drift.
* RK4 shows fourth-order convergence and gives the smallest short-term errors of the three methods.
* Velocity Verlet is second-order accurate and has better bounded energy behavior over long integrations of conservative orbits.
* The Newtonian simulations agree with the analytical Keplerian solution, giving a direct reference for testing the numerical methods.
* The relativistic correction produces periapsis precession, with the leading-order term consistent with $6\pi/L^2$ and the next-order term consistent with $45\pi/L^4$.
* Relativistic effects become stronger at smaller orbital radii.
* Circular orbits become increasingly sensitive to radial perturbations as $r_c$ approaches the ISCO at $r_c=6$.

## Future Improvements

* Study how relativistic periapsis precession changes with eccentricity and angular momentum.
* Investigate zoom-whirl behavior for orbits close to the stability boundary.
* Analyze the Schwarzschild effective potential and connect it to the orbital behavior seen in the simulations.
* Investigate relativistic gravitational scattering and the transition between scattering and capture.
* Extend the model to null geodesics and study gravitational light deflection.
* Study the separatrix between bound, zoom-whirl, scattering, and plunging trajectories.

## Conclusion

This project started with a comparison of numerical methods for Newtonian orbital motion and then extended the same simulation framework to relativistic effects.

The Newtonian simulations use the analytical Keplerian solution as a benchmark. Comparing Euler's Method, RK4, and Velocity Verlet shows the difference between short-term accuracy and long-term stability. RK4 gives the smallest short-term errors, while Velocity Verlet has better bounded energy behavior over long integrations.

I then add the relativistic correction and measure the resulting periapsis precession. The numerical results agree with the predicted leading-order and next-to-leading-order terms, giving a direct way to compare the simulation with the analytical expansion.

Finally, the ISCO analysis shows how circular orbits become increasingly sensitive to radial perturbations as the orbital radius approaches $r_c=6$.

Overall, the project combines numerical methods with orbital mechanics and relativity, using simulations to test both numerical behavior and physical predictions.

## Project Structure

```text
Orbital-Dynamics-and-Relativistic-Precession/
├── docs/
├── images/
├── analysis.py
├── main.py
├── physics.py
├── requirements.txt
├── solvers.py
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
