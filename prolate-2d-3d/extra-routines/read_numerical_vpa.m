% Read data from prolate spheroidal functions into SymPy classes
% Function
%   1. reads file from provided filename in string format
%   2. splits string into substrings with "\n" delimiter
%   3. each substring being a number is converted to SymPy class via 'vpa'
%   4. converted values are stacked in a vertical array

function vpa_array = read_numerical_to_vpa(filename, vpa_digits)

  vpa_digits_orig = digits();

  digits(vpa_digits);
  str_cell_array = strsplit(fileread(filename), "\n")(1:end-1);
  vpa_array = vpa(zeros(length(str_cell_array), 1));

  for i = 1:length(str_cell_array)
    vpa_array(i) = vpa(char(str_cell_array(i)));
  endfor

  digits(vpa_digits_orig);
endfunction
