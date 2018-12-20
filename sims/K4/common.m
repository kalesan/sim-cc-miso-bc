addpath('../..');
addpath('../../cvx');
cvx_setup

SNR = 10 : 5 : 40;

ITERATIONS = 100;
REALIZATIONS = 500;

K = 4;
L = K-1;
