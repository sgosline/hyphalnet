"""
This module contains various statistics needed to compare hyphae specifically
how the communities related to each other and the trees they are built from
"""
import hyphalnet.hypha as hyp
import pandas as pd

def compute_all_distances(hyp_dict):
    """
    Computes distances between all elements of hyphae in list
    Parameters
    ----------
    hyp_dict: Dictionary of hyphae objects

    Returns
    ----------
    Pandas data frame with the following keys
    hyp1: name of hypha 1
    hyp2: name of hypha 2 (can be the same)
    net1: name of network/community 1
    net2: name of network/community 2
    net1_type: type of network for net1
    net2_type: type of network for net2
    distance: distance between two elements, at this point it's jaccard
    """
    df_list = []
    for key1, hyp1 in hyp_dict.items():
        print("Computing distances from", key1)
        #compute forest distances
        within_dist = hyp1.distVals
        within_dist['hyp1'] = key1
        within_dist['hyp2'] = key1
        df_list.append(within_dist)
        for key2, hyp2 in hyp_dict.items():
            if key1 != key2:
                comm_net = hyp1.intra_distance(hyp2)
                comm_net['hyp1'] = key2
                comm_net['hyp2'] = key1
                df_list.append(pd.DataFrame(comm_net))
                #inter_distance is NOT symmetric
                comm_net = hyp2.intra_distance(hyp1)
                comm_net['hyp1'] = key1
                comm_net['hyp2'] = key2
                df_list.append(pd.DataFrame(comm_net))
    dist_df = pd.concat(df_list)
    return dist_df

def compute_all_nmi(hyp_dict, gfile):
    """
    Compute the community assignment between all communities
    """
    print("Computing NMI")
    full_community = hyp.communityFromGraphFile(gfile)
    df_list = []
    for key1, hyp1 in hyp_dict.items():
        net_dist = {'NMI':hyp.computeCommunityNMI(full_community,hyp1.communities)}
        net_dist['hyp1'] = key1
        net_dist['hyp2'] = 'Full Network'
        df_list.append(net_dist)
        for key2,hyp2 in hyp_dict.items():
            if key1 != key2:
                net_dist={'NMI': hyp.computeCommunityNMI(hyp1.communities,hyp2.communities)}
                net_dist['hyp1']=key1
                net_dist['hyp2']=key2
                df_list.append(net_dist)
    return pd.DataFrame(df_list)
