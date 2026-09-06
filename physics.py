import numpy as np


def orbit_derivative(position_and_velocity, t, consts):
    """
    Compute the Newtonian equations of motion for an orbiting body.

    Parameters
    ----------
    position_and_velocity : array_like, shape (4,)
        State vector [x, y, vx, vy].
    t : float
        Current time. Included for compatibility with the solver interface.
    consts : tuple
        Constants passed by the solver. Not used for the normalized
        Newtonian equations implemented here.

    Returns
    -------
    numpy.ndarray, shape (4,)
        Derivative of the state vector [vx, vy, ax, ay].
    """
    x, y, vx, vy = position_and_velocity
    r = np.hypot(x, y)

    ax = -x / r**3
    ay = -y / r**3

    return np.array([vx, vy, ax, ay])


def relativistic_orbit_derivative(position_and_velocity, t, consts):
    """
    Compute the equations of motion including the relativistic correction.

    The effective potential is

        V(r) = -1/r - L²/r³,

    where L is the initial angular momentum. The resulting acceleration is

        a = -r_vec/r³ - 3L² r_vec/r⁵.

    Parameters
    ----------
    position_and_velocity : array_like, shape (4,)
        State vector [x, y, vx, vy].
    t : float
        Current time. Included for compatibility with the solver interface.
    consts : tuple
        Tuple containing the angular momentum L as its first element.

    Returns
    -------
    numpy.ndarray, shape (4,)
        Derivative of the state vector [vx, vy, ax, ay].
    """
    x, y, vx, vy = position_and_velocity
    L = consts[0]

    r = np.hypot(x, y)

    ax = -x / r**3 - 3 * L**2 * x / r**5
    ay = -y / r**3 - 3 * L**2 * y / r**5

    return np.array([vx, vy, ax, ay])


def orbit_acceleration(position, consts):
    """
    Compute the Newtonian gravitational acceleration at a given position.

    Parameters
    ----------
    position : array_like, shape (2,)
        Position vector [x, y].
    consts : tuple
        Constants passed by the solver. Not used for the normalized
        Newtonian equations implemented here.

    Returns
    -------
    numpy.ndarray, shape (2,)
        Gravitational acceleration vector.
    """
    r = np.linalg.norm(position)
    return -position / r**3
