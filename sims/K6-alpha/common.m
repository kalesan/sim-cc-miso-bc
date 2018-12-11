% Common parameters 

addpath('../..');
addpath('../../cvx');

cvx_setup

SNR = -10 : 5 : 40;

ITERATIONS = 100;
realizations = 500;
REALIZATIONS = realizations;

K = 6;   % Number of users
L = K-1; % Number of antennas
S = 2;   % CC size (2 users per message)
