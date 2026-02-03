import pandapower.topology as top
import networkx as nx


# Module-level cache for branch index lookup
_branch_idx_cache = {}


def get_dict_busdir_to_branchidx(net):
    """
    Create mapping from (from_bus, to_bus) to line index.

    Returns unidirectional mapping based on line definition.
    """
    bus_pair_to_line = dict()
    for idx in net.line.index:
        from_bus, to_bus = net.line.from_bus.loc[idx], net.line.to_bus.loc[idx]
        bus_pair_to_line[(from_bus, to_bus)] = idx
    return bus_pair_to_line


def get_branch_variables(net, sending_bus, receiving_bus):
    """
    Get branch variables for power flow analysis.

    Parameters
    ----------
    net : pandapowerNet
        Solved pandapower network
    sending_bus : int
        Sending end bus index
    receiving_bus : int
        Receiving end bus index

    Returns
    -------
    dict
        Dictionary containing branch variables:
        - rp_xq: r*p_out + x*q_out (voltage drop term)
        - z_sq: r^2 + x^2 (impedance squared)
        - z_sq_loss: (r^2 + x^2) * power_loss_ratio
        - r, x: resistance and reactance in pu
        - v_send_sq, v_recv_sq: squared voltage magnitudes at sending/receiving ends
        - p_out, q_out: active/reactive power output in pu
        - power_loss_ratio: (p_out^2 + q_out^2) / v_send_sq
        - z_sq_s_sq: (r^2 + x^2) * (p_out^2 + q_out^2)
        - s_sq: p_out^2 + q_out^2 (apparent power squared)
    """
    global _branch_idx_cache

    vn_kv = float(net.bus.vn_kv.loc[0])
    Zbase = (vn_kv**2) / net.sn_mva

    # Build cache if needed
    net_id = id(net)
    if net_id not in _branch_idx_cache:
        _branch_idx_cache[net_id] = get_dict_busdir_to_branchidx(net)
    bus_pair_to_line = _branch_idx_cache[net_id]

    p_out, q_out, line_idx = None, None, None
    if (sending_bus, receiving_bus) in bus_pair_to_line:
        line_idx = bus_pair_to_line[(sending_bus, receiving_bus)]
        p_out = net.res_line.p_from_mw[line_idx] / net.sn_mva
        q_out = net.res_line.q_from_mvar[line_idx] / net.sn_mva

    elif (receiving_bus, sending_bus) in bus_pair_to_line:
        line_idx = bus_pair_to_line[(receiving_bus, sending_bus)]
        p_out = net.res_line.p_to_mw[line_idx] / net.sn_mva
        q_out = net.res_line.q_to_mvar[line_idx] / net.sn_mva

    else:
        raise Exception(f'There is no such branch {sending_bus}-{receiving_bus}')

    r = net.line.loc[line_idx, 'length_km'] * net.line.loc[line_idx, 'r_ohm_per_km'] / Zbase
    x = net.line.loc[line_idx, 'length_km'] * net.line.loc[line_idx, 'x_ohm_per_km'] / Zbase
    v_send_sq = net.res_bus.loc[sending_bus, 'vm_pu']**2
    v_recv_sq = net.res_bus.loc[receiving_bus, 'vm_pu']**2
    s_sq = p_out**2 + q_out**2
    z_sq = r**2 + x**2
    power_loss_ratio = s_sq / v_send_sq
    rp_xq = r * p_out + x * q_out
    z_sq_loss = z_sq * power_loss_ratio
    z_sq_s_sq = z_sq * s_sq

    return {
        'rp_xq': rp_xq, 'z_sq_loss': z_sq_loss, 'r': r, 'x': x,
        'v_send_sq': v_send_sq, 'v_recv_sq': v_recv_sq,
        'p_out': p_out, 'q_out': q_out,
        'power_loss_ratio': power_loss_ratio, 'z_sq_s_sq': z_sq_s_sq,
        's_sq': s_sq, 'z_sq': z_sq
    }


def get_leaf_buses(net, include_trafos=True, respect_switches=True):
    """Get leaf (terminal) buses in the network, excluding slack bus."""
    G = top.create_nxgraph(
        net,
        include_lines=True,
        include_trafos=include_trafos,
        include_switches=True,
        respect_switches=respect_switches,
        multi=False
    )

    slack = set(net.ext_grid.bus.values)
    leaf_buses = [n for n, deg in G.degree() if deg == 1 and n not in slack]
    return leaf_buses


def path_bus1_to_bus2(net, bus_from, bus_to):
    """Find shortest path between two buses."""
    G = top.create_nxgraph(net)
    path = nx.shortest_path(G, source=bus_from, target=bus_to)
    return path
