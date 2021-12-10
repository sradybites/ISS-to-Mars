function [iCp] = peri2helio(w,I,O)
 dcms = DCMs();
 helio2peri = @(Omega,inc,omega) dcms{3}(omega) * dcms{1}(inc) * dcms{3}(Omega);
 iCp = helio2peri(O,I,w).';
end

