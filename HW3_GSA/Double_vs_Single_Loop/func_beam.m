function y = func_beam(X)
P = X(:, 1);
E = X(:, 2);
I = X(:, 3);

    y = P./E./I;

end