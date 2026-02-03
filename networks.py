import pickle

def star_network():
    with open("networks/starnet.pkl", "rb") as f:
        net = pickle.load(f)
    return net

def twobus_net(v_slack_pu=1.0, r_ohm=0.2, x_ohm=0.4, vn_kv=20.0):
    with open("networks/twobus.pkl", "rb") as f:
        net = pickle.load(f)
    return net

def ieee123():
    with open("networks/ieee123.pkl", "rb") as f:
        net = pickle.load(f)
    return net
