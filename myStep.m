function u = myStep(x,A,h,n)
    u = x;
    u(u < 0) = 0;
    u(u >= h*n) = A;
    u(u >= 0 & u < h*n) = A*u(u >= 0 & u < h*n)/(h*n);
end