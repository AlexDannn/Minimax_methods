%{
Differential correction method.
1) Accepts:
a) Approximate function, given discretely:
• x is a column vector containing a set of points on the approximation interval [-1, 1];
• y is a row vector containing the values ​​of the function at the points of the column vector x;
The number of points with which a function is described affects the number of linear inequalities when solving a linear programming problem. With a large number of such points (500 or more), the method starts to take a lot of time. During testing, it was found that 200 points will be enough to achieve high accuracy.
b) m and n are integers denoting the degrees of the numerator and denominator of the approximating function.

2) Returns:
a) Max_err is the maximum value of the error function on the approximation interval;
b) The graph of the error function and the graph of the approaching and approximating function;
c) a_opt, b_opt are the coefficients of the polynomials of the numerator and denominator of the found approximating function;
d) iter_n - the number of iterations performed by the method.
%}

function [res_max, Errors] = Diff_corr(y, x, m, n)
    Errors = [];
    mmErrors = [];
    if (m >= n)
    TpowersN = []
    TpowersD = []
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
    iter_n = 1;
    
    cvx_begin
        cvx_quiet(true);
        variable a(m+1);
        variable b(n+1);
        subject to
           abs(TpowersN*a - y.*(TpowersD*b)) <= (TpowersD*b);
           abs(b) <= ones(n+1,1);
           TpowersD*b >= 0;
    cvx_end
    r_val = [];
    for i = 1:length(x)
        r_val = [r_val; r_res(x(i), a, b)];
    end
    
    err_func = r_val - y;
    
    [ymax,imax,ymin,imin] = extrema(err_func);
    if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
        mmErrors = [mmErrors;  max(ymax) - min(ymax)];
    else
        mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
    end
    newY = cat(1,ymax,ymin);
    newY = abs(newY);
    
    Errors = [Errors; max(err_func)];

    u=max(newY); l=0; 
    bisection_tol=1e-15; 
    while u-l>= bisection_tol & iter_n < 15
        gamma=(l+u)/2;
        iter_n = iter_n + 1;
        gamma;
        cvx_begin 
            cvx_quiet(true);
            variable a(m+1);
            variable b(n+1);
            subject to
                abs(TpowersN*a-y.*(TpowersD*b)) <= gamma*TpowersD*b;
                abs(b) <= ones(n+1,1);
                TpowersD*b >= 0;
        cvx_end

        if strcmp(cvx_status,'Solved')
            u=gamma;
            a_opt=a;
            b_opt=b;
            objval_opt=gamma;
            y_fit=TpowersN*a_opt./(TpowersD*b_opt);
            
            Errors = [Errors; max(y_fit-y)];
            err_func = y_fit - y;
            
            [ymax,imax,ymin,imin] = extrema(err_func);
            if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
                mmErrors = [mmErrors;  max(ymax) - min(ymax)];
            else
                mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
            end
            
        else
            a_opt = a;
            b_opt = b;
            y_fit=TpowersN*a_opt./(TpowersD*b_opt);
            Errors = [Errors; max(y_fit-y)];
            l=gamma;
            err_func = y_fit - y;
            
            [ymax,imax,ymin,imin] = extrema(err_func);
            if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
                mmErrors = [mmErrors;  max(ymax) - min(ymax)];
            else
                mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
            end
        end


    
    
    end
        figure(2*iter_n+1);
    plot(x, y_fit-y, 'Linewidth',2);
    set(gca,'FontSize',18)
    legend('Error function'),grid on
    disp([vpa(mmErrors)]);
    figure(2*iter_n);
    semilogy(mmErrors, 'Linewidth',2);
    set(gca,'FontSize',18)
    legend('Minimax Error'),grid on;
    end

    if (m < n)
    TpowersN = []
    TpowersD = []
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
    iter_n = 1
    
    cvx_begin
        cvx_quiet(true);
        variable a(m+1);
        variable b(n+1);
        subject to
           abs(TpowersN*a - y.*(TpowersD*b)) <= (TpowersD*b);
           abs(b) <= ones(n+1,1);
           TpowersD*b >= 0;
        cvx_end
        
        
    r_val = []
    for i = 1:length(x)
        r_val = [r_val; r_res(x(i), a, b)];
    end
    
    err_func = r_val - y    
    [ymax,imax,ymin,imin] = extrema(err_func);
    newY = cat(1,ymax,ymin);
    newY = abs(newY);
    disp([min(newY)]);
    disp([max(newY)]);
    
    
    u=max(newY); l=0;
    bisection_tol=1e-15;
    while u-l>= bisection_tol & iter_n < 15
        gamma=(l+u)/2;
        iter_n = iter_n + 1
        gamma;
        cvx_begin
            cvx_quiet(true);
            variable a(m+1);
            variable b(n+1);
            subject to
                abs(TpowersN*a-y.*(TpowersD*b)) <= gamma*TpowersD*b;
                abs(b) <= ones(n+1,1);
                TpowersD*b >= 0;
        cvx_end
    
        if strcmp(cvx_status,'Solved')
            u=gamma;
            a_opt=a;
            b_opt=b;
            objval_opt=gamma;
        else
            a_opt = a;
            b_opt = b;
            l=gamma;
        end
    y_fit=TpowersN*a_opt./(TpowersD*b_opt);

    figure(2*iter_n);
    plot(x,y,'b', x,y_fit,'r');
    grid on
    legend('F(x)','R(x)')

    figure(2*iter_n+1);
    plot(x, y_fit-y);
    grid on
    legend('Error func')
    end
    end
    a_opt = a;
    b_opt = b;

    res_max = max(y_fit-y)
end



