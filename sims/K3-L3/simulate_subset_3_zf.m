common;

cc(precoder_zf(), 'data/sim-subset-3-zf.mat', L, K, S, SNR, ...
   'REALIZATIONS', REALIZATIONS, 'name', 'CC-ZF (power loading)');
