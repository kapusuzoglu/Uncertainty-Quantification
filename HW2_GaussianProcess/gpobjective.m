function f = gpobjective(x,xi,Y)

k = eye(length(Y));
for i=1:length(Y)
    for j=1:length(Y)
        if (i<j)
            k(i,j)= (-0.5*(((xi(i,1)-xi(j,1))/x(1))^2 ...
            + ((xi(i,2)-xi(j,2))/x(2))^2 +  ((xi(i,3)-xi(j,3))/x(3))^2));
        elseif (j<i)
            k(i,j) = k(j,i);
        end
    end
end
k = x(4)^2 * exp(k);

% D = size(x,1);

%   % 3. Create the contribution due to squared exponential ARD.    
%       KMN = pdist2(xi(:,1)/x(1), xi(:,1)/x(1)).^2;
%       for r = 2:D-1
%           KMN = KMN + pdist2(xi(:,r)/x(r), xi(:,r)/x(r)).^2;        
%       end
%       KMN = (x(4)^2)*exp(-0.5*KMN)
      
% f = -0.5*Y'*( (k+eye(length(Y))*sigma_e^2)\Y) - 0.5*log(det(k+eye(length(Y))*...
%     sigma_e^2)) - length(Y)/2*log(2*pi);
f = -0.5*Y'*(k\Y) - 0.5*log(det(k)) - length(Y)/2*log(2*pi);
% argmax_x f = argmin_x f
f = -f;

end