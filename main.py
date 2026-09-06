import numpy as np
import matplotlib.pyplot as plt

from physics import (
    orbit_derivative,
    orbit_acceleration,
    relativistic_orbit_derivative,
)
from solvers import euler_solver, rk4_solver, velocity_verlet_solver
from analysis import (
    get_analytical_orbit,
    get_analytical_period,
    get_precession_amplitude,
)


initial_position = (3, 0)


# ============================================================
# GRAPH 1 — Numerical convergence
# ============================================================

timesteps = np.linspace(0.002, 0.1, 20)
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# Initial conditions chosen to compare eccentric and near-circular orbits.
for col, Vy0 in enumerate([0.3, 0.55]):
    initial_velocity = (0, Vy0)
    initial_state = initial_position + initial_velocity

    energy = Vy0**2 / 2 - 1 / 3
    L = 3 * Vy0
    eccentricity = np.sqrt(1 + 2 * energy * L**2 / 1**2)
    orbit_type = "Eccentric" if Vy0 == 0.3 else "Near-Circular"
    title = fr"{orbit_type} Orbit: $e={eccentricity:.3f}$"

    errors = {"Euler": [], "RK4": [], "Velocity Verlet": []}

    for dt in timesteps:
        times = np.arange(0, 30 + dt, dt)

        euler = euler_solver(orbit_derivative, initial_state, times)
        rk4 = rk4_solver(orbit_derivative, initial_state, times)
        positions, _ = velocity_verlet_solver(
            orbit_acceleration,
            initial_position,
            initial_velocity,
            times
        )

        x_exact, y_exact = get_analytical_orbit(
            initial_position,
            initial_velocity,
            [times[-1]],
            1,
        )

        exact_position = np.array([x_exact[0], y_exact[0]])

        errors["Euler"].append(
            np.linalg.norm(euler[-1, :2] - exact_position)
        )
        errors["RK4"].append(
            np.linalg.norm(rk4[-1, :2] - exact_position)
        )
        errors["Velocity Verlet"].append(
            np.linalg.norm(positions[-1] - exact_position)
        )

    slopes = {
        method: np.polyfit(np.log(timesteps), np.log(error), 1)[0]
        for method, error in errors.items()
    }

    # Compare convergence on linear and log-log scales.
    axes[0, col].plot(timesteps, errors["RK4"], label="RK4")
    axes[0, col].plot(
        timesteps,
        errors["Velocity Verlet"],
        label="Velocity Verlet",
    )
    axes[0, col].set(
        xlabel="Timestep Size",
        ylabel="Position Error",
        title=title,
    )
    axes[0, col].legend()
    axes[0, col].grid()

    for method, error in errors.items():
        axes[1, col].loglog(
            timesteps,
            error,
            label=f"{method} (slope = {slopes[method]:.2f})",
        )

    axes[1, col].set(
        xlabel="Timestep",
        ylabel="Position Error",
        title=title,
    )
    axes[1, col].legend()
    axes[1, col].grid()

fig.suptitle(
    "Numerical Convergence for Newtonian Orbital Integration",
    fontsize=16,
)
fig.tight_layout()
fig.subplots_adjust(top=0.92)
plt.savefig("numerical_convergence.png")
plt.show()


# ============================================================
# GRAPH 2 — Long-term energy conservation
# ============================================================

fig, axes = plt.subplots(3, 2, figsize=(12, 10))

for col, Vy0 in enumerate([0.3, 0.55]):
    initial_velocity = (0, Vy0)
    initial_state = initial_position + initial_velocity

    energy0 = Vy0**2 / 2 - 1 / 3
    L = 3 * Vy0
    eccentricity = np.sqrt(1 + 2 * energy0 * L**2 / 1**2)
    orbit_type = "Eccentric" if Vy0 == 0.3 else "Near-Circular"
    title = fr"{orbit_type} Orbit: $e={eccentricity:.3f}$"

    period = get_analytical_period(initial_position, initial_velocity, 1)
    dt = 0.1

    # Euler and RK4 are compared over 500 periods.
    times = np.arange(0, 500 * period + dt, dt)

    for row, solver in enumerate([euler_solver, rk4_solver]):
        states = solver(orbit_derivative, initial_state, times)

        r = np.linalg.norm(states[:, :2], axis=1)
        speed = np.linalg.norm(states[:, 2:], axis=1)
        energy = speed**2 / 2 - 1 / r
        relative_error = (energy - energy[0]) / abs(energy[0]) * 100

        axes[row, col].plot(times, relative_error)
        axes[row, col].set(
            xlabel="Time",
            ylabel="Relative Energy Error (%)",
            title=title,
        )
        axes[row, col].grid()

    if Vy0 == 0.55:
        axes[1, col].yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x / 1e-5:.0f} × 10⁻⁵")
        )

    # Velocity Verlet only needs a shorter interval for comparison.
    times = np.arange(0, 10 * period + dt, dt)
    positions, velocities = velocity_verlet_solver(
        orbit_acceleration,
        initial_position,
        initial_velocity,
        times
    )

    r = np.linalg.norm(positions, axis=1)
    speed = np.linalg.norm(velocities, axis=1)
    energy = speed**2 / 2 - 1 / r
    relative_error = (energy - energy[0]) / abs(energy[0]) * 100

    axes[2, col].plot(times, relative_error)
    axes[2, col].set(
        xlabel="Time",
        ylabel="Relative Energy Error (%)",
        title=title,
    )
    axes[2, col].grid()

for row, label in enumerate(["Euler", "RK4", "Velocity Verlet"]):
    axes[row, 0].text(
        -0.25,
        0.5,
        label,
        transform=axes[row, 0].transAxes,
        rotation=90,
        va="center",
        ha="center",
        fontsize=14,
    )

fig.suptitle(
    "Long-Term Energy Conservation of Numerical Orbital Integrators",
    fontsize=16,
)
fig.tight_layout()
fig.subplots_adjust(top=0.93, left=0.10)
plt.savefig("energy_conservation.png")
plt.show()


# ============================================================
# GRAPH 3 — Relativistic orbital dynamics
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, Vy0, total_periods, dr, dt_factor, title in [
    (
        axes[0], 0.01, 30, 0.1, 0.01,
        "Orbital Dynamics in the Weak-Field Limit",
    ),
    (
        axes[1], 0.15, 3, 0.1, 0.01,
        "Relativistic Precession at Intermediate Radius",
    ),
    (
        axes[2], 0.6, None, 0.001, 0.001,
        "Stable and Unstable Radial Perturbations Near the ISCO",
    ),
]:
    x0 = 1 / Vy0**2 + 3
    unperturbed_state = np.array((x0, 0, 0, Vy0))

    period = 2 * np.pi / Vy0**3
    L = x0 * Vy0

    perturbations = [(dr, 0, 0, 0), (-dr, 0, 0, 0)]

    # Near the ISCO, the two perturbations require different runtimes.
    periods = [10, 3.79] if total_periods is None else [total_periods] * 2

    for perturbation, orbit_periods, label in zip(
        perturbations,
        periods,
        ["Outward perturbation $+\\epsilon$",
         "Inward perturbation $-\\epsilon$"],
    ):
        initial_state = unperturbed_state * (1 + np.array(perturbation))
        times = np.arange(
            0,
            orbit_periods * period + dt_factor * period,
            dt_factor * period,
        )

        states = rk4_solver(
            relativistic_orbit_derivative,
            initial_state,
            times,
            (L,),
        )

        ax.plot(states[:, 0], states[:, 1], label=label)

    ax.plot(0, 0, "ko", ms=8)
    ax.set(
        xlabel="x-position",
        ylabel="y-position",
        title=title,
    )
    ax.grid()
    ax.axis("equal")

fig.suptitle(
    "Relativistic Orbital Dynamics and Stability",
    fontsize=16,
)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="lower center",
    ncol=2,
    bbox_to_anchor=(0.5, -0.02),
)

fig.tight_layout()
fig.subplots_adjust(top=0.85, bottom=0.15)
plt.savefig("relativistic_orbits.png")
plt.show()


# ============================================================
# GRAPH 4 — Relativistic precession coefficients
# ============================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig3, ax3 = plt.subplots(figsize=(7, 5))

epsilons = np.linspace(0.01, 0.1, 7)
velocities = np.linspace(0.03, 0.45, 21) ** 2
epsilon_characteristic = 0.01

for epsilon in epsilons:
    L_values = []
    precession_values = []

    for Vy0 in velocities:
        x0 = (1 + epsilon) / Vy0**2 + 3
        initial_state = (x0, 0, 0, Vy0)

        period = 2 * np.pi / Vy0**3
        L = x0 * Vy0
        times = np.arange(
            0,
            10 * period + 0.0005 * period,
            0.0005 * period,
        )

        states = rk4_solver(
            relativistic_orbit_derivative,
            initial_state,
            times,
            (L,),
        )

        L_values.append(L)
        precession_values.append(
            get_precession_amplitude(states, times)
        )

    L_values = np.array(L_values)
    precession_values = np.array(precession_values)

    x = 1 / L_values**2
    y = precession_values * L_values**2

    order = np.argsort(x)
    x, y = x[order], y[order]

    intercepts = []
    slopes = []

    for n_fit in range(3, 10):
        slope, intercept = np.polyfit(x[:n_fit], y[:n_fit], 1)
        slopes.append(slope)
        intercepts.append(intercept)

    label = fr"$\epsilon = {epsilon:.2f}$"

    ax1.plot(range(3, 10), intercepts, marker="o", label=label)
    ax2.plot(range(3, 10), slopes, marker="o", label=label)

    if np.isclose(epsilon, epsilon_characteristic):
        characteristic_L = L_values
        characteristic_precession = precession_values

ax1.axhline(6 * np.pi, linestyle="--", label="Theoretical: $6\\pi$")
ax1.set(
    xlabel="Number of points used in fit",
    ylabel="Fitted leading-order coefficient",
    title="Leading-Order Precession Coefficient",
    ylim=(18.7, 19.0),
)
ax1.legend()

ax2.axhline(45 * np.pi, linestyle="--", label="Theoretical: $45\\pi$")
ax2.set(
    xlabel="Number of points used in fit",
    ylabel="Fitted next-to-leading order coefficient",
    title="Next-to-Leading-Order Precession Coefficient",
    ylim=(139, 144),
)
ax2.legend()

fig.suptitle(
    "Numerical Extraction of Relativistic Precession Coefficients",
    fontsize=16,
)
fig.tight_layout()
fig.subplots_adjust(top=0.85)

x = 1 / characteristic_L**2
ax3.plot(
    x,
    characteristic_precession,
    "o",
    label=r"Numerical, $\epsilon=0.01$",
)
ax3.plot(
    x,
    6 * np.pi / characteristic_L**2,
    "--",
    label=r"$6\pi/L^2$",
)
ax3.plot(
    x,
    6 * np.pi / characteristic_L**2
    + 45 * np.pi / characteristic_L**4,
    label=r"$6\pi/L^2 + 45\pi/L^4$",
)

ax3.set(
    xlabel=r"$1/L^2$",
    ylabel=r"Precession $\Delta\phi$ (rad/orbit)",
    title=r"Numerical vs. Theoretical Precession ($\epsilon=0.01$)",
)
ax3.legend()
fig3.tight_layout()
plt.savefig("relativistic_precession.png")
plt.show()


# ============================================================
# GRAPH 5 — Critical perturbation near the ISCO
# ============================================================

Vy0_values = np.linspace(0.53, 0.63, 21)
epsilon_values = np.logspace(
    np.log10(0.00005),
    np.log10(0.005),
    20,
)

for periods in [10, 15, 20, 25]:
    critical_epsilons = []

    for Vy0 in Vy0_values:
        x0 = 1 / Vy0**2 + 3
        unperturbed_state = np.array((x0, 0, 0, Vy0))
        period = 2 * np.pi / Vy0**3

        epsilon_crit = None

        # Find the smallest perturbation that produces a >10% radial deviation.
        for epsilon in epsilon_values:
            perturbations = np.array([
                (epsilon, 0, 0, 0),
                (-epsilon, 0, 0, 0),
            ])

            times = np.arange(
                0,
                periods * period + 0.01 * period,
                0.01 * period,
            )

            for initial_state in unperturbed_state * (1 + perturbations):
                x_init, y_init, Vx_init, Vy_init = initial_state
                L = x_init * Vy_init - y_init * Vx_init

                states = rk4_solver(
                    relativistic_orbit_derivative,
                    initial_state,
                    times,
                    (L,),
                )

                radius = np.linalg.norm(states[:, :2], axis=1)
                radial_error = (radius - x0) / x0

                if np.any(np.abs(radial_error) > 0.1):
                    epsilon_crit = epsilon
                    break

            if epsilon_crit is not None:
                break

        if epsilon_crit is not None:
            critical_epsilons.append(epsilon_crit)

    radius_values = 1 / Vy0_values**2 + 3

    plt.plot(
        radius_values,
        critical_epsilons,
        label=fr"$t_{{max}} = {periods}T$",
    )

plt.axvline(
    6,
    linestyle="--",
    label=r"Theoretical ISCO: $r_c=6$",
)

plt.xlabel(r"Circular orbit radius $r_c$")
plt.ylabel(r"Critical fractional perturbation $\epsilon_{\mathrm{crit}}$")
plt.title("Critical Perturbation for Orbits Near the ISCO")
plt.legend()
plt.grid()
plt.savefig("isco_critical_perturbation.png")
plt.show()