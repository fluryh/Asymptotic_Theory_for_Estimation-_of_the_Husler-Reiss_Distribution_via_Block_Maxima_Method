function out = rhr(n,l)
out = rand(n,2);
out(:,1) = -log(-log(out(:,1)));
for j = 1:n
    f = @(y) condDist(y,out(j,2),l,out(j,1));
    out(j,2) = fminunc(f,2);
end
end