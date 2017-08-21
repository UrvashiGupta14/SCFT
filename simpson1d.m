function s = simpson1d( y, dx )
n = length(y) - 1;
s = dx/3*(y(1)+4*sum(y(2:2:n))+2*sum(y(3:2:n-1))+y(n+1));
end