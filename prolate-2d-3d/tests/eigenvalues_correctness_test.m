% NOTE : to be run from 'spheroidal folder'
% before to run PSWFs must be re-generated via legacy code of G. Sabinin
%


eigvals = zeros(1, 30);

for ind = 0:29
  ind
  pswf_struct = pro_open_everything('saved', 10, 0, ind);
  pswf = pswf_struct.S1.S1;
  fc_pswf = TransformFcSimpson1d(pswf, linspace(-1,1, 1025), 10, linspace(-1,1,1025));
  rand_ind = randi(1025);
  eigvals(ind + 1) = fc_pswf(rand_ind) ./ pswf(rand_ind);

  printf("pswf index %d :  legacy %s : check %s : diff %f\n", ind, ...
                                                  num2str(pswf_struct.coefficients.mu), ...
                                                  num2str(eigvals(ind+1)), ...
                                                  abs(pswf_struct.coefficients.mu - eigvals(ind+1)));
endfor
