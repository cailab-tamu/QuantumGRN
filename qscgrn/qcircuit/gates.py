import numpy as np


__all__ = ["ry_gate", "der_ry_gate", "cry_gate", "der_cry_gate"]


def ry_gate(theta):
    """
    Computes the matrix transformation for the `Ry` gate with an
    angle parameter of `theta`.
    `T = Ry(theta)`.
    Parameters
    ----------
    theta : float or int
        The rotation angle in radians. In q-sphere, the rotation angle
        around the y-axis.
    Returns
    -------
    T : ndarray
    Raises
    ------
    TypeError
        If the input is not an float or integer.
    """
    if not isinstance(theta, (float, int)):
        raise TypeError("Ry gate operation unsupported for {type}"
                        .format(type=type(theta)))

    half_theta = theta / 2
    sin = np.sin(half_theta)
    cos = np.cos(half_theta)
    return np.array([[cos, -sin],
                     [sin, cos]])


def der_ry_gate(theta):
    """
    Computes the matrix transformation for the derivative of `Ry` gate
    with an angle parameter of `theta`.
    `dT = dRy(theta) dtheta `.
    Parameters
    ----------
    theta : float or int
        The rotation angle in radians. In q-sphere, the rotation angle
        around the y-axis.
    Returns
    -------
    dT : ndarray
    Raises
    ------
    TypeError
        If tthe input is not an float or integer.
    """
    if not isinstance(theta, (float, int)):
        raise TypeError("dRy gate operation unsupported for {type}"
                        .format(type=type(theta)))

    half_theta = theta / 2
    sin = np.sin(half_theta)
    cos = np.cos(half_theta)
    return 0.5 * np.array([[-sin, -cos],
                           [cos, -sin]])


def cry_gate(theta, nqubits, control, target):
    """
    Computes the matrix transformation for the `c-Ry` gate with an
    angle parameter of `theta`.
    `T = c-Ry(theta)`.

    The quantum circuit notation follows the qiskit convention. In
    other words: a quantum state is written as `|x_3x_2x_1x_0>`.

    q0 --------------.-----------------
                     |
    q1 --------------+-----------------
                     |
    q2 --------------+-----------------
                     |
    q3 --------------Ry----------------

    Parameters
    ----------
    theta : float or int
        The rotation angle in radians. In q-sphere, the rotation angle
        around the y-axis.
    nqubits : int
        The number of qubits in the quantum circuit. Similar
        to the number of genes for modelling.
    control : int
        The control qubit, which activates or not the `Ry` gate in the
        target qubit.
    target : int
        The target qubit, where the `Ry` gate will be applied.
    Returns
    -------
    T : ndarray
    Raises
    ------
    TypeError
        If the input theta is not an float or integer.
    """
    if not isinstance(theta, (float, int)):
        raise TypeError("c-Ry gate operation unsupported for {type}"
                        .format(type=type(theta)))

    zero_matrix = np.array([[1., 0.], [0., 0.]])
    one_matrix = np.array([[0., 0.], [0., 1.]])

    RR, II = None, None

    for i in range(nqubits):
        if i == control:
            RR = np.kron(one_matrix, RR) if RR is not None \
                else one_matrix
            II = np.kron(zero_matrix, II) if II is not None \
                else zero_matrix

        elif i == target:
            RR = np.kron(ry_gate(theta), RR) if RR is not None \
                else ry_gate(theta)
            II = np.kron(np.identity(2), II) if II is not None \
                else np.identity(2)

        else:
            RR = np.kron(np.identity(2), RR) if RR is not None \
                else np.identity(2)
            II = np.kron(np.identity(2), II) if II is not None \
                else np.identity(2)

    return RR + II


def der_cry_gate(theta, nqubits, control, target):
    """
    Computes the matrix transformation for the derivative of `c-Ry`
    gate with an angle parameter of `theta`.
    `dT = dc-Ry(theta) dtheta`.
    Parameters
    ----------
    theta : float or int
        The rotation angle in radians. In q-sphere, the rotation angle
        around the y-axis.
    nqubits : int
        The number of qubits in the quantum circuit. Similar
        to the number of genes for modelling.
    control : int
        The control qubit, which activates or not the `Ry` gate in the
        target qubit.
    target : int
        The target qubit, where the `Ry` gate will be applied.
    Returns
    -------
    dT : ndarray
    Raises
    ------
    TypeError
        If the input theta is not an float or integer.
    """
    if not isinstance(theta, (float, int)):
        raise TypeError("c-Ry gate operation unsupported for {type}"
                        .format(type=type(theta)))
    zero_matrix = np.array([[1., 0.], [0., 0.]])
    one_matrix = np.array([[0., 0.], [0., 1.]])

    RR = None

    for i in range(nqubits):
        if i == control:
            RR = one_matrix if i == 0 else np.kron(one_matrix, RR)

        elif i == target:
            RR = der_ry_gate(theta) if i == 0 \
                else np.kron(der_ry_gate(theta), RR)

        else:
            RR = np.identity(2) if i == 0 \
                else np.kron(np.identity(2), RR)

    return RR

