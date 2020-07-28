F1 = @(z) sin(z);
F2 = @(z) sin(2.5*cos(z));
F3 = @(z) sin(z).*cos(z);
F4 = @(z) cos(4*sin(z));
F5 = @(z) exp(sin(z));

F6 = @(z) exp(z);
F7 = @(z) sinh(z);
F8 = @(z) atan(z);
F9 = @(z) z.^6;
F10 = @(z) z.^7 + 3*z.^5 + 2*z -1;
F11 = @(z) log(z+3);
F12 = @(x) 1./(1+x.^2);
F16 = @(z) sin(z)./(z.^2+1);
F18 = @(z) max(sin(20*z), exp(z-1));
F19 = @(z) tanh(z + 0.5) - tanh(z - 0.5);

F20 = @(z) log(1.001 + z);
F21 = @(z) (100*pi*(z.^2-0.36))./(sinh(100*pi*(z.^2-0.36)));
F22 = @(z) abs(z).*sqrt(abs(z));
F23 = @(z) abs(z);
F24 = @(z) tanh(50*z);
F25 = @(z) sin(20*z);
F26 = @(x) (x-(x.^3)/factorial(3)+ (x.^5)/factorial(5))./(x.^2+1);

        
X = linspace(-1, 1, 1000);
tic
 [r,pol,res,ZER, ZJ, FJ, WJ, ERRVEC, wt, mmErrors] = aaa(F24,X,'degree',22)
times = vpa(toc)


figure(1)
plot(X, F24(X),X, r(X), 'Linewidth', 2), set(gca,'FontSize',18), grid on
legend('r(x)', 'f(x)'), grid on
 
figure(2)
plot(X, F24(X) - r(X)), set(gca,'FontSize',18)
legend('Error function'),grid on


% n_iter = []; 
% iter_err = [];
% for i = 1:length(ERRVEC)
%   n_iter = [n_iter; i-1];
%   iter_err = [iter_err; ERRVEC(i)];
% end


vpa(ERRVEC)

  
