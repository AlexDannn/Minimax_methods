        F1 = @(z) sin(z);
        F2 = @(z) sin(2.5*cos(z));
        F3 = @(z) sin(z).*cos(z);
        F4 = @(z) cos(sin(z));
        F5 = @(z) exp(sin(z));


        F6 = @(z) exp(z);
        F7 = @(z) sinh(z);
        F8 = @(z) atan(z);
        F9 = @(z) z.^6;
        F10 = @(z) z.^7 + 3*z.^5 + 2*z -1;
        F11 = @(z) log(z+3);

        F12 = @(z) 1./(z.^2+1);
        F13 = @(z) (z.^5 + z.^3 + z)./(z.^6 + z.^2 + 3);
        F14 = @(z) 1./(z.^2 + z + 1);
        F15 = @(z) sin(z)./z;
        F16 = @(z) sin(z)./(z.^2+1);
        F17 = @(z) (z.^2 + 2)./(z.^2 + 1);


        F18 = @(z) max(sin(20*z), exp(z-1));
        F19 = @(z) tanh(z + 0.5) - tanh(z - 0.5);
        F20 = @(z) log(1.001 + z);
        F21 = @(z) (100*pi*(z.^2-0.36))./(sinh(100*pi*(z.^2-0.36)));
        F22 = @(z) abs(z).*sqrt(abs(z));
        F23 = @(z) abs(z);
        F24 = @(z) tanh(50*z);
        F25 = @(z) sin(20*z);
        

        m = 11;
        n = 11;
        x = linspace(-1, 1, 271)';
        y = abs(x).*sqrt(abs(x))
        
        
       tic
       [res_max, Errors] = Diff_corr(y, x, m, n);
       Exec_time = toc

