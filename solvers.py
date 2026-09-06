import numpy as np


def euler_solver(derivative_function, initial_state, times, consts=()):
    """
    Solve a system of ordinary differential equations using Euler's method.

    Parameters
    ----------
    derivative_function : callable
        Function with signature ``f(state, time, consts)`` that returns
        the derivative of the state.
    initial_state : array_like
        Initial state vector.
    times : array_like
        Times at which the solution is evaluated.
    consts : tuple, optional
        Additional constants passed to the derivative function.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(len(times), len(initial_state))`` containing
        the state at each time.
    """
    states = np.zeros((len(times), len(initial_state)))
    states[0] = initial_state

    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]
        derivative = derivative_function(states[i], times[i], consts)
        states[i + 1] = states[i] + derivative * dt

    return states


def rk4_solver(derivative_function, initial_state, times, consts=()):
    """
    Solve a system of ordinary differential equations using fourth-order
    Runge-Kutta integration (RK4).

    Parameters
    ----------
    derivative_function : callable
        Function with signature ``f(state, time, consts)`` that returns
        the derivative of the state.
    initial_state : array_like
        Initial state vector.
    times : array_like
        Times at which the solution is evaluated.
    consts : tuple, optional
        Additional constants passed to the derivative function.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(len(times), len(initial_state))`` containing
        the state at each time.
    """
    states = np.zeros((len(times), len(initial_state)))
    states[0] = initial_state

    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]
        t = times[i]
        state = states[i]

        k1 = derivative_function(state, t, consts)

        k2 = derivative_function(
            state + 0.5 * dt * k1,
            t + 0.5 * dt,
            consts,
        )

        k3 = derivative_function(
            state + 0.5 * dt * k2,
            t + 0.5 * dt,
            consts,
        )

        k4 = derivative_function(
            state + dt * k3,
            t + dt,
            consts,
        )

        states[i + 1] = state + dt / 6 * (
            k1 + 2 * k2 + 2 * k3 + k4
        )

    return states


def velocity_verlet_solver(
    acceleration_function,
    initial_position,
    initial_velocity,
    times,
    consts=(),
):
    """
    Solve a second-order system using the Velocity Verlet method.

    This method is particularly useful for orbital dynamics because it
    has good long-term energy and phase-space behavior for conservative
    systems.

    Parameters
    ----------
    acceleration_function : callable
        Function with signature ``a(position, consts)`` that returns
        the acceleration at a given position.
    initial_position : array_like
        Initial position vector.
    initial_velocity : array_like
        Initial velocity vector.
    times : array_like
        Times at which the position and velocity are evaluated.
    consts : tuple, optional
        Additional constants passed to the acceleration function.

    Returns
    -------
    positions : numpy.ndarray
        Array of shape ``(len(times), len(initial_position))`` containing
        the position at each time.
    velocities : numpy.ndarray
        Array of shape ``(len(times), len(initial_velocity))`` containing
        the velocity at each time.
    """
    positions = np.zeros((len(times), len(initial_position)))
    velocities = np.zeros((len(times), len(initial_velocity)))

    positions[0] = initial_position
    velocities[0] = initial_velocity

    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]

        acceleration = acceleration_function(positions[i], consts)

        positions[i + 1] = (
            positions[i]
            + velocities[i] * dt
            + 0.5 * acceleration * dt**2
        )

        new_acceleration = acceleration_function(
            positions[i + 1],
            consts,
        )

        velocities[i + 1] = (
            velocities[i]
            + 0.5 * (acceleration + new_acceleration) * dt
        )

    return positions, velocities