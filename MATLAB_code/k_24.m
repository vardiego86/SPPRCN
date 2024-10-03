function rate = k_24(sigma)

kdaj = 2.136e-4;      % Different from values used for RCN model out of multiscaling
kMdaj = 2.136e-4;

rate = -kdaj*sigma + kMdaj;