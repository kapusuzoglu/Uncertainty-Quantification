function gS = get_generalized_indices(S, D, timevec)

gS(1) = 0;
Nt = length(timevec);
for n = 2 : Nt
   ti = timevec(1:n);
   Di = D(1:n);
   Si = S(1:n);
   
   intD = trapz(ti, Di);

   w = Di ./ intD;

   gS(n) = trapz(ti, w .* Si);
end
