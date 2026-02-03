import numpy as np

from powerflow import complex_bus_voltage_pu, bus_injection_current_pu, get_Zbus_reduced_pu
from topology import get_branch_variables, get_leaf_buses, path_bus1_to_bus2


def compute_L_index(net, *, gen_buses=None, load_buses=None):
    """
    Compute Kessel & Glavitsch (1986) L-index for voltage stability.

    Parameters
    ----------
    net : pandapowerNet
        Solved pandapower network
    gen_buses : list, optional
        Generator bus indices. If None, uses ext_grid and gen buses.
    load_buses : list, optional
        Load bus indices. If None, uses all non-generator buses.

    Returns
    -------
    L_by_bus : dict
        L-index for each load bus
    L_max : float
        Maximum L-index value
    crit_bus : int
        Bus with maximum L-index
    """
    Ybus = net._ppc["internal"]["Ybus"]
    Y = Ybus.tocsc()

    pp2ppc = net._pd2ppc_lookups["bus"]

    if gen_buses is None:
        G_pp = set(int(b) for b in net.ext_grid.bus.values)
        if len(net.gen):
            G_pp |= set(int(b) for b in net.gen.bus.values)
    else:
        G_pp = set(int(b) for b in gen_buses)

    all_pp = set(int(b) for b in net.bus.index)
    if load_buses is None:
        L_pp = list(sorted(all_pp - G_pp))
    else:
        L_pp = list(int(b) for b in load_buses)

    G_pp = list(sorted(G_pp))

    if len(G_pp) == 0 or len(L_pp) == 0:
        raise ValueError("Generator bus set or load bus set is empty.")

    G = [int(pp2ppc[b]) for b in G_pp]
    L = [int(pp2ppc[b]) for b in L_pp]

    Y_LL = Y[L, :][:, L].toarray()
    Y_LG = Y[L, :][:, G].toarray()

    F = -np.linalg.solve(Y_LL, Y_LG)

    V_bus = complex_bus_voltage_pu(net)
    n = Y.shape[0]
    V = np.zeros(n, dtype=complex)
    for b in net.bus.index:
        V[int(pp2ppc[int(b)])] = V_bus[int(b)]

    Vg = V[G]

    L_by_bus = {}
    for row, b_pp in enumerate(L_pp):
        i = int(pp2ppc[b_pp])
        Vi = V[i]
        if abs(Vi) < 1e-12:
            Li = np.inf
        else:
            Li = abs(1.0 - np.sum(F[row, :] * (Vg / Vi)))
        L_by_bus[int(b_pp)] = float(Li)

    crit_bus = max(L_by_bus, key=lambda b: L_by_bus[b])
    L_max = L_by_bus[crit_bus]
    return L_by_bus, float(L_max), int(crit_bus)


def compute_margins(net, slack_bus):
    """
    Compute multiple voltage stability margins.

    Parameters
    ----------
    net : pandapowerNet
        Solved pandapower network
    slack_bus : int
        Slack bus index

    Returns
    -------
    inj_based_margin : dict
        Wang/Cui/Wang injection-based margin for each bus
    L_by_bus : dict
        Kessel & Glavitsch L-index for each load bus
    single_branch_det : dict
        Single-branch determinant margin for each line
    multiple_branch_deri : dict
        Path-accumulated margin based on derivative
    """
    # --- Wang/Cui/Wang margin ---
    V = complex_bus_voltage_pu(net)
    Iinj = bus_injection_current_pu(net, V)
    Zred, keep = get_Zbus_reduced_pu(net, slack_bus_idx=slack_bus)
    Ired = Iinj[keep]
    inj_based_margin = {}
    for a_red, h in enumerate(keep):
        rhs = 0.0
        for b_red, i_bus in enumerate(keep):
            rhs += np.abs(Zred[a_red, b_red] * Ired[b_red])
        inj_based_margin[h] = np.abs(V[h]) - rhs

    # --- Kessel & Glavitsch (1986) L-index ---
    L_by_bus, Lmax, Lcrit = compute_L_index(net)

    # --- Single branch level margin ---
    single_branch_det = {}
    for line_idx in net.line.index:
        branch_vars = get_branch_variables(net, net.line.from_bus[line_idx], net.line.to_bus[line_idx])
        r, x = branch_vars['r'], branch_vars['x']
        p_out, q_out = branch_vars['p_out'], branch_vars['q_out']
        v_send_sq = branch_vars['v_send_sq']
        single_branch_det[int(line_idx)] = (v_send_sq - 2*(r*p_out + x*q_out))**2

        branch_vars = get_branch_variables(net, net.line.to_bus[line_idx], net.line.from_bus[line_idx])
        r, x = branch_vars['r'], branch_vars['x']
        p_out, q_out = branch_vars['p_out'], branch_vars['q_out']
        v_send_sq = branch_vars['v_send_sq']
        single_branch_det[len(net.line.index) + int(line_idx)] = (v_send_sq - 2*(r*p_out + x*q_out))**2

    buses = list(set(net.bus.index))

    # --- Multiple branch level margin based on derivative ---
    multiple_branch_deri = dict()
    multiple_branch_deri[net.ext_grid.bus[0]] = 999

    for bus in buses:
        multiple_branch_deri[(bus, slack_bus)] = accumulate_determinant(net, bus, slack_bus)

    return inj_based_margin, L_by_bus, single_branch_det, multiple_branch_deri


def accumulate_determinant(net, bus_from, bus_to, terminate_at_slack=True):
    """
    Compute path-accumulated voltage stability determinant.

    Parameters
    ----------
    net : pandapowerNet
        Solved pandapower network
    bus_from : int
        Starting bus
    bus_to : int
        Ending bus
    terminate_at_slack : bool
        If True, terminate path at slack bus

    Returns
    -------
    float
        Determinant value (999 if path length < 2)
    """
    path = path_bus1_to_bus2(net, bus_from, bus_to)

    if len(path) < 2:
        return 999
    if (net.ext_grid.bus[0] in path) and (bus_from != net.ext_grid.bus[0]):
        bus_to = net.ext_grid.bus[0]

    path = path_bus1_to_bus2(net, bus_from, bus_to)
    sum_rp_xq = 0
    sum_power_term = 0
    voltage_sensitivity = dict()
    voltage_ratio = dict()

    v_from_sq, v_to_sq = net.res_bus.vm_pu[bus_from]**2, net.res_bus.vm_pu[bus_to]**2

    for i in range(len(path)-1):
        branch_vars = get_branch_variables(net, path[i], path[i+1])
        rp_xq = branch_vars['rp_xq']
        z_sq = branch_vars['z_sq']
        s_sq = branch_vars['s_sq']
        v_send_sq = branch_vars['v_send_sq']

        if i == 0:
            voltage_sensitivity[path[i]] = 1

        sum_rp_xq += rp_xq
        sum_power_term += z_sq * s_sq * v_from_sq / v_send_sq
        voltage_ratio[(path[i], path[i+1])] = (1 - z_sq * s_sq / (v_send_sq**2))
        voltage_sensitivity[path[i+1]] = voltage_sensitivity[path[i]] * voltage_ratio[(path[i], path[i+1])]

    determinant = (v_to_sq + 2*sum_rp_xq)**2 - 4*sum_power_term
    return determinant
