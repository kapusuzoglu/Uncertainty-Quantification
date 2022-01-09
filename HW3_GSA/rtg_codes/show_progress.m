function show_proress(ratio)
if mod(ratio, 5)<1e-13
   if ratio == 5
      fprintf('[');
   end
   fprintf('%g%% ', ratio);
   if ratio == 100
      fprintf(']\n');
   end
end
