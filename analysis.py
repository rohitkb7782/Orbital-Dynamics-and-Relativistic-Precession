import numpy as np
from scipy.optimize import fsolve


def get_eccentric_anomaly(t, a, e, GM=1, C=0):
    """
    Solve Kepler's equation for the eccentric anomaly at time t.

    Kepler's equation is

        E - e sin(E) = n t + C,

    where n is the mean motion.

    Parameters
    ----------
    t : float
        Time.
    a : float
        Semi-major axis of the orbit.
    e : float
        Orbital eccentricity.
    GM : float, optional
        Gravitational parameter. Defaults to 1.
    C : float, optional
        Phase constant determined by the initial conditions.

    Returns
    -------
    float
        Eccentric anomaly E at time t.
    """
    n = np.sqrt(GM / a**3)
    mean_anomaly = n * t + C

    def kepler_equation(E):
        return E - e * np.sin(E) - mean_anomaly

    return fsolve(kepler_equation, mean_anomaly)[0]


def get_analytical_orbit(initial_position, initial_velocity, times, GM):
    """
    Calculate the analytical Newtonian orbit using Keplerian motion.

    The initial conditions are used to determine the orbital energy,
    angular momentum, semi-major axis, eccentricity, and phase.

    Parameters
    ----------
    initial_position : array_like, shape (2,)
        Initial position vector [x, y].
    initial_velocity : array_like, shape (2,)
        Initial velocity vector [vx, vy].
    times : array_like
        Times at which the analytical position is evaluated.
    GM : float
        Gravitational parameter.

    Returns
    -------
    x_analytical : numpy.ndarray
        Analytical x-position at each time.
    y_analytical : numpy.ndarray
        Analytical y-position at each time.
    """
    x0, y0 = initial_position
    Vx0, Vy0 = initial_velocity

    r0 = np.linalg.norm(initial_position)
    speed0 = np.linalg.norm(initial_velocity)

    energy = speed0**2 / 2 - GM / r0
    Lz = x0 * Vy0 - y0 * Vx0

    a = -GM / (2 * energy)
    e = np.sqrt(1 + 2 * energy * Lz**2 / GM**2)

    E0 = np.arctan2(
        y0 / (a * np.sqrt(1 - e**2)),
        e - x0 / a,
    )
    C = E0 - e * np.sin(E0)

    E_values = np.array([
        get_eccentric_anomaly(t, a, e, GM, C)
        for t in times
    ])

    x_analytical = a * (e - np.cos(E_values))
    y_analytical = -a * np.sqrt(1 - e**2) * np.sin(E_values)

    return x_analytical, y_analytical


def get_analytical_period(initial_position, initial_velocity, GM):
    """
    Calculate the orbital period from the initial conditions.

    Parameters
    ----------
    initial_position : array_like, shape (2,)
        Initial position vector [x, y].
    initial_velocity : array_like, shape (2,)
        Initial velocity vector [vx, vy].
    GM : float
        Gravitational parameter.

    Returns
    -------
    float
        Orbital period.
    """
    r0 = np.linalg.norm(initial_position)
    speed0 = np.linalg.norm(initial_velocity)

    energy = speed0**2 / 2 - GM / r0
    a = -GM / (2 * energy)

    return 2 * np.pi * np.sqrt(a**3 / GM)


def get_precession_amplitude(states, times):
    """
    Determine the average periapsis-to-periapsis angular advance.

    Periapsides are identified from changes in the radial velocity
    from negative to positive. Their angular positions are then
    unwrapped before calculating the mean change in angle.

    Parameters
    ----------
    states : numpy.ndarray, shape (N, 4)
        Orbit states [x, y, vx, vy] at each time.
    times : array_like, shape (N,)
        Times corresponding to the states.

    Returns
    -------
    float
        Mean angular separation between successive periapsides in radians.
    """
    x, y = states[:, 0], states[:, 1]
    vx, vy = states[:, 2], states[:, 3]

    r = np.hypot(x, y)
    radial_velocity = (x * vx + y * vy) / r

    periapsis_angles = []

    # Detect periapsis by locating where radial velocity changes from negative to positive.
    for i in range(1, len(times)):
        if radial_velocity[i - 1] < 0 <= radial_velocity[i]:
            # Interpolate the zero crossing of radial velocity.
            fraction = (
                -radial_velocity[i - 1]
                / (radial_velocity[i] - radial_velocity[i - 1])
            )

            x_peri = x[i - 1] + fraction * (x[i] - x[i - 1])
            y_peri = y[i - 1] + fraction * (y[i] - y[i - 1])

            periapsis_angles.append(np.arctan2(y_peri, x_peri))

    periapsis_angles = np.unwrap(periapsis_angles)

    return np.mean(np.diff(periapsis_angles))