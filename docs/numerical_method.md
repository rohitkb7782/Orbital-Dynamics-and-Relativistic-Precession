## Numerical Methods

I compare three numerical integrators: Euler's Method, fourth-order Runge-Kutta (RK4), and Velocity Verlet.

### Euler's Method

Euler's Method was also used in my previous project, [*Projectile Motion with Quadratic Drag*](https://github.com/rohitkb7782/Projectile-Motion-with-Quadratic-Drag), where its first-order convergence was compared with an analytical solution. It advances the state using the derivative at the current timestep:

$$
\mathbf s_{n+1}
=\mathbf s_n+
\mathbf f(\mathbf s_n,t_n)h.
$$

It has global error $O(h).$ Here, I use Euler's Method mainly as a baseline for comparison.

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

RK4 has global error $O(h^4).$

### Velocity Verlet

Velocity Verlet is useful for systems where the acceleration depends on position. It first updates the position using the current velocity and acceleration:

$$
\mathbf r_{n+1}
=\mathbf r_n+
\mathbf v_nh+
\frac12\mathbf a_nh^2.
$$

The acceleration is then recalculated from the new position:

$$
\mathbf a_{n+1}
=\mathbf a(\mathbf r_{n+1}),
$$

and the velocity is updated using the average of the two accelerations:

$$
\mathbf v_{n+1}
=\mathbf v_n+
\frac12
(\mathbf a_n+\mathbf a_{n+1})h.
$$

Velocity Verlet has global error $O(h^2),$ and its time-symmetric structure gives it good long-term energy behavior for conservative systems.

### Convergence and Energy Conservation

I test numerical convergence by comparing the numerical position with the Keplerian reference solution using the relative position error:

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

The expected scaling is $\epsilon_r\propto h^p,$ where $p$ is the order of the method.

For long-term behavior, I track the relative energy error:

$$
\frac{E(t)-E(0)}{|E(0)|}\times100\%.
$$

This lets me compare the short-term accuracy of the methods with how well they conserve energy over longer simulations.
