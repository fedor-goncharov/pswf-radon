
sqnorms = zeros(1, 30);

for ind = 0:29
  ind
  pswf_struct = pro_open_everything('saved', 10, 0, ind);
  pswf = pswf_struct.S1.S1;
  sqnorms(ind+1) = IntegrateSimpson1D(pswf .* pswf, 2.0/1024);

  printf("squared norm : pswf index %d :  theoretical %1.1f : check %f : diff %f\n", ind, ...
                                                  1.0, ...
                                                  sqnorms(ind+1), ...
                                                  abs(1.0-sqnorms(ind+1)));
endfor

% conclusion : norms are nearly 1 on [-1,1]
