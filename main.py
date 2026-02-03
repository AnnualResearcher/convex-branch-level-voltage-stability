import numpy as np
import pandapower as pp
import matplotlib.pyplot as plt

from networks import ieee123, twobus_net, star_network
from plotting import set_plot_style
from powerflow import get_svd
from metrics import compute_margins


def sweep_any_net(net, lam_values, network_name="star"):
    """
    Sweep load multiplier and compute voltage stability metrics.

    Parameters
    ----------
    net : pandapowerNet
        Pandapower network to analyze
    lam_values : array-like
        Load multiplier values to sweep
    network_name : str
        Name for output file labeling
    """
    slack_bus = net.ext_grid.bus[0]

    base_load_p = net.load.p_mw.values.copy()
    base_load_q = net.load.q_mvar.values.copy()
    base_sgen_p = net.sgen.p_mw.values.copy()
    base_sgen_q = net.sgen.q_mvar.values.copy()
    min_singular_values = []
    results = []

    for load_mult in lam_values:
        net.load.p_mw = load_mult * base_load_p
        net.load.q_mvar = load_mult * base_load_q
        net.sgen.p_mw = load_mult * base_sgen_p
        net.sgen.q_mvar = load_mult * base_sgen_q
        try:
            pp.runpp(net, algorithm="nr", init="flat", calculate_voltage_angles=True,
                     max_iteration=50, tolerance_mva=1e-8)
            min_svd = get_svd(net)
            min_singular_values.append(min_svd)
        except Exception:
            results.append((load_mult, False, np.nan, np.nan, None, None, None, None, None, None))
            continue

        inj_margin, l_index, single_branch, path_accumulated = compute_margins(net, slack_bus=slack_bus)

        inj_margin_min = float(np.min(list(inj_margin.values())))
        l_index_max = float(np.max(list(l_index.values())))
        single_branch_min = float(np.min(list(single_branch.values())))
        path_accumulated_min = float(np.min(list(path_accumulated.values())))

        inj_margin_crit_bus = int(min(inj_margin, key=lambda b: inj_margin[b]))
        l_index_crit_bus = int(min(l_index, key=lambda b: l_index[b]))
        single_branch_crit_line = int(min(single_branch, key=lambda line: single_branch[line]))
        path_accumulated_crit = min(path_accumulated, key=lambda b: path_accumulated[b])

        results.append((load_mult, True, inj_margin_min, l_index_max, single_branch_min, path_accumulated_min,
                        inj_margin_crit_bus, l_index_crit_bus, single_branch_crit_line, path_accumulated_crit))

    print("load_mult, pf_converged, inj_margin_min, l_index_max, single_branch_min, path_accum_min, crit_buses_lines")
    for row in results:
        print(row)

    try:
        set_plot_style(column="one", font_size=8)

        load_multipliers = np.array([row[0] for row in results], dtype=float)
        pf_converged = np.array([row[1] for row in results], dtype=bool)

        inj_margin_arr = np.array([row[2] for row in results], dtype=float)
        l_index_arr = np.array([row[3] for row in results], dtype=float)
        single_branch_arr = np.array([row[4] for row in results], dtype=float)
        path_accumulated_arr = np.array([row[5] for row in results], dtype=float)

        x_valid = load_multipliers[pf_converged]
        inj_margin_valid = inj_margin_arr[pf_converged]
        l_index_valid = l_index_arr[pf_converged]
        single_branch_valid = single_branch_arr[pf_converged]
        path_accumulated_valid = path_accumulated_arr[pf_converged]

        fig, ax = plt.subplots()

        ax.plot(x_valid, l_index_valid, color="grey", linewidth=1.0, linestyle="-", label="L-index [1]")
        ax.plot(x_valid, inj_margin_valid, color="black", linewidth=1.0, linestyle="-", label=r"$|V|-\sum |ZI|$ [2]")

        if network_name == 'two_bus':
            single_branch_linestyle = (0, (5, 5))
            path_accumulated_linestyle = (5, (5, 5))
        else:
            single_branch_linestyle = '-'
            path_accumulated_linestyle = '-'

        ax.plot(x_valid, single_branch_valid, color="#0E21A0", linewidth=1.0,
                linestyle=single_branch_linestyle, label="Single-branch")
        ax.plot(x_valid, path_accumulated_valid, color="#d62728", linewidth=1.0,
                linestyle=path_accumulated_linestyle, label="Path-accumulated")

        ax.set_xlabel(r"Load multiplier $\lambda$")
        ax.set_ylabel("Metric")

        ax.grid(True, which="major", linewidth=0.6, alpha=0.35)
        ax.grid(False, which="minor")
        ax.minorticks_off()

        ax.tick_params(direction="in", top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

        ax.legend(frameon=False)

        fig.savefig(f"figures/metric_{network_name}.svg", format='svg')
        plt.show()

    except ImportError:
        pass


if __name__ == "__main__":
    import argparse

    NETWORK_CONFIGS = {
        'twobus': {
            'loader': twobus_net,
            'lam_max': 309.01698,
            'name': 'two_bus',
            'description': 'Two-bus network'
        },
        'ieee123': {
            'loader': ieee123,
            'lam_max': 3.04914,
            'name': 'ieee',
            'description': 'IEEE 123-bus test feeder'
        },
        'star': {
            'loader': star_network,
            'lam_max': 26.21818,
            'name': 'star',
            'description': 'Star network'
        }
    }

    parser = argparse.ArgumentParser(
        description='Branch-level voltage stability analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python main.py twobus      Run analysis on two-bus network
  python main.py ieee123     Run analysis on IEEE 123-bus test feeder
  python main.py star        Run analysis on star network
  python main.py all         Run analysis on all networks
        '''
    )
    parser.add_argument(
        'network',
        choices=['twobus', 'ieee123', 'star', 'all'],
        help='Network to analyze: twobus, ieee123, star, or all'
    )
    parser.add_argument(
        '--points', '-n',
        type=int,
        default=200,
        help='Number of load multiplier points (default: 200)'
    )

    args = parser.parse_args()

    if args.network == 'all':
        networks_to_run = ['twobus', 'ieee123', 'star']
    else:
        networks_to_run = [args.network]

    for network_key in networks_to_run:
        config = NETWORK_CONFIGS[network_key]
        print(f"\n{'='*60}")
        print(f"Running: {config['description']}")
        print(f"{'='*60}")

        net = config['loader']()
        lam_values = np.linspace(1, config['lam_max'], args.points)
        sweep_any_net(net, lam_values, network_name=config['name'])
