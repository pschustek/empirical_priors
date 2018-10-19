% Fig. 3
mu = [0.54 0.28 0.41 0.37];
N = [7 5 12 3]';
nH = [3 1 5 2]';

[ confHeads, bEv, pbEv, confBlockHeads ] = opt_inf.all( nH, N, [1 1 1 1]', 14, 9 );