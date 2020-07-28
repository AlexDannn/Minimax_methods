F1 = @(x) sin(x);
F2 = @(x) sin(2.5*cos(x));
F3 = @(x) sin(x).*cos(x);
F4 = @(x) cos(sin(x));
F5 = @(x) exp(sin(x));

F6 = @(x) exp(x);
F7 = @(z) sinh(x);
F8 = @(x) atan(x);
F9 = @(x) x.^6;
F10 = @(x) x.^7 + 3*x.^5 + 2*x -1;
F11 = @(x) log(x+3);

F12 = @(x) 1./(x.^2+1);
F13 = @(x) (x.^5 + x.^3 + x)./(x.^6 + x.^2 + 3);
F14 = @(x) 1./(x.^2 + x + 1);
F15 = @(x) sin(x)./x;
F16 = @(x) sin(x)./(x.^2+1);
F17 = @(x) (x.^2 + 2)./(x.^2 + 1);

F18 = @(x) max(sin(20*x), exp(x-1));
F19 = @(x) tanh(x + 0.5) - tanh(x - 0.5);
F20 = @(x) log(1.001 + x);
F21 = @(x) (100*pi*(x.^2-0.36))./(sinh(100*pi*(x.^2-0.36)));
F22 = @(x) abs(x).*sqrt(abs(x));
F23 = @(x) abs(x);
F24 = @(x) tanh(50*x);


%{
DCR2 method.
1) Accepts:
a) Approximate function, given discretely:
• x is a column vector containing a set of points on the approximation interval [-1, 1];
• y is a row vector containing the values of the function at the points of the column vector x;
The number of points is formed based on the entered degrees m and n (m + n + 2);
b) m and n are integers denoting the degrees of the numerator and denominator of the approximating function.

2) Returns:
a) Max_err is the maximum value of the error function on the approximation interval;
b) The graph of the error function and the graph of the approaching and approximating function;
c) a, b are the coefficients of the polynomials of the numerator and denominator of the found approximating function;
d) iter_n - the number of iterations performed by the method.
%}

m = 17; n = 3;
f = @(x) tanh(50*x);

x = linspace(-1,1,m+n+2)';
z = linspace(-1,1,2001)';
y = f(x);
y_z = f(z);

errVec = []
mmErrors = []; 

TpowersD = zeros(length(x), n+1);
TpowersN = zeros(length(x), m+1);
    for i = 1:length(x)
        for j=1:m+1
            TpowersN(i,j) = x(i)^(j-1);
        end
    end
   
     for i = 1:length(x)
         for j=1:n+1
             TpowersD(i,j) = x(i)^(j-1);
        end
     end
     
    cvx_begin
        cvx_quiet(true);
        variable a(m+1);
        variable b(n+1);
        subject to
            abs(TpowersN*a - y.*(TpowersD*b)) <= (TpowersD*b);
            abs(b) <= ones(n+1,1);
            TpowersD*b >= 0;
    cvx_end
    
    aPrev = a;
    bPrev = b;
    r_val = [];
    
    for i = 1:length(z)
        r_val = [r_val; r_res(z(i), aPrev, bPrev)];
    end
    
    err_func = r_val - y_z;
    
    [ymax,imax,ymin,imin] = extrema(err_func);
    if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
        mmErrors = [mmErrors;  max(ymax) - min(ymax)];
    else
        mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
    end

    [ymax,imax,ymin,imin] = extrema(err_func);


stopDCR = 0;
X_sets = []';
min_max_rev = [];
error_abs = 1;
err1 = max(abs(err_func));
tic
iter_n = 0;

while error_abs > 10e-15 & stopDCR < 35
     iter_n = iter_n + 1;
     errVec = [errVec; err1];
     newX = [];
     yNew = [];
     
     newX = cat(1,z(imax),z(imin));
     newX = sort(newX);
     X_prev = newX;
     X_sets = [X_sets; newX];
     yNew = [];
     for i=1:length(newX)
        yNew = [yNew; f(newX(i))];
     end
     
     TpowersD = zeros(length(newX), n+1);
     TpowersN = zeros(length(newX), m+1);
 
     for i = 1:length(newX)
             for j=1:m+1
                 TpowersN(i,j) = newX(i)^(j-1);
             end
         end
 
          for i = 1:length(newX)
              for j=1:n+1
                  TpowersD(i,j) = newX(i)^(j-1);
             end
          end

         cvx_begin 
         cvx_quiet(true);
         variable a(m+1);
         variable b(n+1);
         subject to
             abs(TpowersN*a - yNew.*(TpowersD*b)) <= (err1*(TpowersD*b));
             abs(b) <= ones(n+1,1);
             TpowersD*b >= 0;

         cvx_end
     r_valNew = [];
    
     for i = 1:length(z)
        r_valNew = [r_valNew; r_res(z(i), a, b)];
     end
     err_func = [];
     err_func = r_valNew - y_z;
     imax = [];
     imin = [];
     ymax = [];
     ymin = [];
     
    [ymax,imax,ymin,imin] = extrema(err_func);
    if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
        mmErrors = [mmErrors;  max(ymax) - min(ymax)];
    else
        mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
    end
    
     checkX = cat(1,z(imax),z(imin));
     checkX = sort(checkX);
     checkY = [];
     
     for i = 1:length(checkX)
         checkY = [checkY f(checkX(i)) - r_res(checkX(i), a, b)];
     end
     
     errPrev = err1;
     err1 = max(abs(err_func));

     stopDCR = stopDCR + 1;
     min_min = abs(max(ymin));
     max_max = abs(max(ymax));
     min_max_rev = [min_max_rev; abs(min_min)/(max_max) - 1]; 
     error_abs = abs(min_min - max_max);

end 
  figure
  plot(z, err_func)
  legend('err function')
  grid on
  
  
Ex_time = toc;  
disp([max(r_valNew - y_z)]);
iter_n
a, b
