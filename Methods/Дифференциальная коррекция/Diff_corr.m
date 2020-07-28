%{
����� ���������������� ���������. 
1)	���������:
    a)	 ������������ �������, �������� ���������:
        x ? ������-�������, ���������� ����� ����� �� ��������� ������������� [-1, 1]; 
        y ? ������-������, ���������� �������� ������� � ������ �������-������� x;
���������� �����, � ������� ������� ����������� �������, ������ �� ���������� �������� ���������� ��� ������� ������ ��������� ����������������. ��� ������� ���������� ����� ����� (500 � �����) ������ ������ �������� �������� ����� �������. � ���� ������������ ���� ��������, ��� 200 ����� ����� ���������� ��� ���������� ������� ��������.  
    �)	 m � n ? ����� �����, ������������ ������� ��������� � ����������� ������������ �������.

2)	����������:
    �)	Max_err ? ������������ �������� ������� ������ �� ��������� �������������;
    �)	������ ������� ������ � ������ ������������ � ������������ �������;
    �)	a_opt, b_opt ? ������������ ��������� ��������� � ����������� ��������� ������������ �������;
    �)	iter_n ? ���������� ��������, ����������� �������.
%}

function [res_max, Errors] = Diff_corr(y, x, m, n)
    Errors = [];
    mmErrors = [];
    % ������, ����� ������� ��������� ������, ��� ������� �����������.
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
    
    % ����� �������� �����������.
    cvx_begin
        cvx_quiet(true);
        variable a(m+1);
        variable b(n+1);
        subject to
           abs(TpowersN*a - y.*(TpowersD*b)) <= (TpowersD*b);
           abs(b) <= ones(n+1,1);
           TpowersD*b >= 0;
    cvx_end
    % ���������� ������������ �������.
    r_val = [];
    for i = 1:length(x)
        r_val = [r_val; r_res(x(i), a, b)];
    end
    
    err_func = r_val - y;
    
    % ���������� ����������� ������.
    [ymax,imax,ymin,imin] = extrema(err_func);
    if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
        mmErrors = [mmErrors;  max(ymax) - min(ymax)];
    else
        mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
    end
    newY = cat(1,ymax,ymin);
    newY = abs(newY);
    
    Errors = [Errors; max(err_func)];
    % �������� ���� ��������� � �������������� ������ ��������.
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
            % ���������� ������������ �������. 
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
    % ������, ����� ������� ����������� ������, ��� ������� ���������.
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
    % ������ ������������ � ������������ �������.
    figure(2*iter_n);
    plot(x,y,'b', x,y_fit,'r');
    grid on
    legend('F(x)','R(x)')
    % ������ ������� ������.
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



