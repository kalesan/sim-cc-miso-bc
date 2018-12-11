common;

cc(precoder_zf(), 'data/sim-subset-6-zf.mat', L, K, S, SNR, ...
   'REALIZATIONS', realizations, 'name', 'CC-ZF');
