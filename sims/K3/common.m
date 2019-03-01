addpath('../..');
addpath('../../cvx');
cvx_setup

SNR = 0 : 5 : 30;

ITERATIONS = 100;
REALIZATIONS = 500;

K = 3;
L = K-1;

S = 2; % Two users per message
