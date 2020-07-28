F1 = @(x) sin(x);
F2 = @(x) sin(2.5*cos(x));
F3 = @(x) sin(x).*cos(x);
F4 = @(x) cos(sin(4*x));
F5 = @(x) exp(sin(x));

F6 = @(x) exp(x);
F7 = @(x) sinh(x);
F8 = @(x) atan(x);
F9 = @(x) x.^6;
F11 = @(x) log(x+3);

F12 = @(x) 1./(x.^2+1);
F13 = @(x) (x.^5 + x.^3 + x)./(x.^6 + x.^2 + 3);
F14 = @(x) 1./(x.^2 + x + 1);
F16 = @(x) sin(x)./(x.^2+1);
F17 = @(x) (x.^2 + 2)./(x.^2 + 1);

F18 = @(x) max(sin(20*x), exp(x-1));
F20 = @(x) log(1.001 + x);
F21 = @(x) (100*pi*(x.^2-0.36))./(sinh(100*pi*(x.^2-0.36)))
F22 = @(x) abs(x).*sqrt(abs(x));
F24 = @(x) tanh(50*x); 
F25 = @(x) sin(20*(x));


%{
Basic IRLS algorithm.
1) Accepts:
a) Approximate function, given discretely:
• x - a list containing a set of points on the approximation interval [-1, 1];
• y - a list containing the values ​​of the function at points in the list x;
The number of points is formed based on the degree of the approximating polynomial N (N + 2);

b) N is an integer denoting the degree of the approximating polynomial (the default range is from 4 to 25).

2) Returns:
a) result - a list of method errors when approximating the approximating function, depending on the degree N;
b) The graph of the error function and the graph of the approaching and approximating function;
c) bestN - the degree of the best approximation function.
%}


LF = 1000; N = 16; 
x = linspace(-1, 1, LF)';


y = (x.^2 + 2)./(x.^2 + 1);
% ---------------------


A = zeros(LF,N);
 for i = 1:LF
     for j = 1:N
         A(i,j) = x(i)^(j-1);
     end
 end
 
 KK = 20;
 p = 2.5;
 [b, iter_n] = Basic_IRLS(A, y, p, KK);
 
 r_x = []
 for i =1:length(x)
     r_x = [r_x; func_sum(x(i),b)];
 end
 

 plot(x, y-r_x), grid on
 

 disp([max(y-r_x)]);
 
 

 function val = func_sum(x, b)
    val = 0;
    for i = 1:length(b)
        val = val + x^(i-1)*b(i);
    end
 end
 
    
