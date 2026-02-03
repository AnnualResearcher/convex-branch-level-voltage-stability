import numpy as np
from numpy.linalg import inv


def complex_bus_voltage_pu(net):
    """Get complex bus voltages in per-unit from solved network."""
    vm = net.res_bus.vm_pu.values
    va = np.deg2rad(net.res_bus.va_degree.values)
    return vm * np.exp(1j * va)


def bus_injection_current_pu(net, V_pu):
    """
    Calculate bus injection currents in per-unit.

    Note: pandapower res_bus p_mw/q_mvar are net injection (gen - load) by sign convention.
    """
    S_pu = (net.res_bus.p_mw.values + 1j * net.res_bus.q_mvar.values) / net.sn_mva
    I_pu = np.conj(S_pu / V_pu)
    return I_pu


def get_Zbus_reduced_pu(net, slack_bus_idx):
    """
    Get reduced Z-bus matrix (slack bus eliminated) in per-unit.

    Note: runpp must have been called before using this function.

    Returns
    -------
    Zred : ndarray
        Reduced impedance matrix
    keep : list
        Mapping from reduced index to original bus index
    """
    Ybus = net._ppc["internal"]["Ybus"]  # scipy sparse (pu)
    Y = Ybus.tocsc()

    n = Y.shape[0]
    keep = [i for i in range(n) if i != slack_bus_idx]
    Yred = Y[keep, :][:, keep].toarray()
    Zred = inv(Yred)
    return Zred, keep


def line_z_pu(net, line_idx):
    """Get series impedance of a line in per-unit."""
    fb = int(net.line.from_bus.loc[line_idx])
    vn_kv = float(net.bus.vn_kv.loc[fb])
    Zbase = (vn_kv**2) / net.sn_mva
    r = float(net.line.r_ohm_per_km.loc[line_idx]) * float(net.line.length_km.loc[line_idx])
    x = float(net.line.x_ohm_per_km.loc[line_idx]) * float(net.line.length_km.loc[line_idx])
    z_ohm = r + 1j * x
    z_pu = z_ohm / Zbase
    return z_pu


def get_svd(net):
    """Get minimum singular value of the Jacobian matrix."""
    J = net._ppc["internal"]["J"]
    a = J.tocsc()
    b = a.toarray()
    return np.linalg.svdvals(b)[-1]
