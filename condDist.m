function val = condDist(y,q,l,x)
    val = abs(exp(exp(-x)*(1-normcdf(l+(y-x)/(2*l))) - exp(-y)*normcdf(l+(x-y)/(2*l)))*normcdf(l+(y-x)/2*l) - q);
end