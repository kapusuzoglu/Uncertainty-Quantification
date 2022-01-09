function output = beam_RBD( input, Mu, Std )

P = norminv(input(:, 1), Mu(1), Std(1));
E = norminv(input(:, 2), Mu(2), Std(2));
I = norminv(input(:, 3), Mu(3), Std(3));

output = func_beam([P, E, I]);


end