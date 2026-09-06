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

For Newtonian gravity,

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

I use these quantities to construct the analytical Keplerian reference solution.

### Keplerian Analytical Solution

The Newtonian orbit can be written using the eccentric anomaly $E_{\rm anom}$, which satisfies Kepler's equation

$$
E_{\rm anom}-e\sin E_{\rm anom}=nt+C,
$$

where $C$ contains the initial conditions and

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

I solve Kepler's equation numerically with `fsolve`. This gives an independent reference solution that I can use to test the numerical integrators.

### Relativistic Binet Equation

A Binet equation is a useful way to describe the shape of an orbit under a central force. It relates the orbit $u=1/r$ to the radial acceleration:

$$
u''+u=-\frac{F(u)}{L^2u^2}.
$$

For the relativistic acceleration,

$$
F(r)=-\frac{1}{r^2}-\frac{3L^2}{r^4}.
$$

Using $u=1/r$,

$$
F(u)=-u^2-3L^2u^4.
$$

Substituting into the Binet equation,

$$
u''+u
=-\frac{-u^2-3L^2u^4}{L^2u^2}
=\frac{1}{L^2}+3u^2.
$$

Therefore, the relativistic Binet equation is

$$
\boxed{u''+u=\frac{1}{L^2}+3u^2}.
$$

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

To study its stability, I introduce a small radial perturbation,

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

When $\omega_r^2>0$,

$$
\delta u\propto\cos(\omega_r\phi),
$$

so the perturbation remains bounded and the orbit is stable.

When $\omega_r^2<0$, the frequency is imaginary and

$$
\delta u\propto e^{\pm\gamma\phi},
$$

so the perturbation grows and the orbit is unstable.

At $\omega_r=0$, the radial restoring force vanishes, giving marginal stability.

Therefore,

$$
r_c>6\Rightarrow\text{stable},
\qquad
r_c=6\Rightarrow\text{marginally stable},
\qquad
r_c<6\Rightarrow\text{unstable}.
$$

This gives the analytical ISCO radius

$$
\boxed{r_{\rm ISCO}=6}.
$$

To initialize the circular orbits used in the simulations, I write the circularity condition in terms of the tangential velocity. The required centripetal acceleration satisfies

$$
\frac{V_t^2}{r}
=\frac1{r^2}+\frac{3L^2}{r^4}.
$$

Using $L=rV_t$ gives

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

I then test the stability of these orbits by applying small radial perturbations,

$$
r_0\rightarrow r_0(1\pm\epsilon),
$$

and measuring the resulting radial deviation.

### Relativistic Precession

The Binet equation

$$
u''+u=\frac1{L^2}+3u^2
$$

is nonlinear because of the $3u^2$ term. In the weak-field limit, $u=1/r$ is small, so the relativistic correction is also small. This makes it possible to treat the relativistic effect as a perturbation of the Newtonian orbit.

I expand the precession in powers of the small parameter $1/L^2$:

$$
\frac{c_1}{L^2}
+
\frac{c_2}{L^4}
+
\mathcal O(L^{-6}),
$$

where $c_1$ is the leading-order coefficient and $c_2$ is the next-to-leading-order coefficient.

For a circular orbit, the radial frequency is

$$
\omega_r^2=1-6u_0,
$$

so one complete radial cycle requires

$$
\Phi=\frac{2\pi}{\omega_r}
=\frac{2\pi}{\sqrt{1-6u_0}}.
$$

In the Newtonian case, the radial cycle closes after one revolution, so $\Phi=2\pi$. The relativistic periapsis advance is therefore

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

For $1/L^2\ll1$, this can be expanded as

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

This gives me a way to test the numerical results. I can extract the coefficients of the $1/L^2$ and $1/L^4$ terms from the simulations and compare them with the theoretical values.

Numerically, I identify periapsides when

$$
v_r=\frac{\mathbf r\cdot\mathbf v}{r}
$$

crosses from negative to positive. The angular separation between successive periapsides gives the measured precession, which I then fit against $1/L^2$.
