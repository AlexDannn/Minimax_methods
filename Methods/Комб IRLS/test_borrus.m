% ������������� �������
F1 = @(x) sin(x);
F2 = @(x) sin(2.5*cos(x));
F3 = @(x) sin(x).*cos(x);
F4 = @(x) cos(sin(4*x));
F5 = @(x) exp(sin(x));

% ���������� ������� 
F6 = @(x) exp(x);
F7 = @(x) sinh(x);
F8 = @(x) atan(x);
F9 = @(x) x.^6;
F11 = @(x) log(x+3);

% ������-������������ �������
F12 = @(x) 1./(x.^2+1);
F13 = @(x) (x.^5 + x.^3 + x)./(x.^6 + x.^2 + 3);
F14 = @(x) 1./(x.^2 + x + 1);
F16 = @(x) sin(x)./(x.^2+1);
F17 = @(x) (x.^2 + 2)./(x.^2 + 1);

% �������������� �������
F18 = @(x) max(sin(20*x), exp(x-1));
F20 = @(x) log(1.001 + x);
F21 = @(x) (100*pi*(x.^2-0.36))./(sinh(100*pi*(x.^2-0.36)));
F22 = @(x) abs(x).*sqrt(abs(x));
F24 = @(x) tanh(50*x); 
F25 = @(x) sin(20*(x));
F26 = @(x) (x-(x.^3)/factorial(3)+ (x.^5)/factorial(5))./(x.^2+1);

%{
��������������� �������� IRLS.
1)	���������:
    �) ������������ �������, �������� ���������:
    	x  ������, ���������� ����� ����� �� ��������� ������������� [-1, 1];
    	y  ������, ���������� �������� ������� � ������ ������ x;
    	p  ����������� �� �������� p;
    	K  �������� �������������;
    ���������� ����� ����������� ������ �� ������� ������������� �������� N (N + 2); 

    �) N  ����� �����, ������������ ������� ������������� �������� (�� ��������� ����� �������� �� 4 �� 25).

2)	����������:
    �)	������ ������� ������ � ������ ������������ � ������������ �������;
    �)	Min_err  ������������ ������ ������� �����������;
    �)  �������� P_k.

%}


iter = 0;
Errors = [];
Exec_times = [];
%  N - ������� ������������ �������. 
for N = 22:22
    mmErrors = [];
    tic
    % ���������� �����, ����������� �� p, �������� ���������� K. 
    LF = 1000; p = 10e5; K = 1.01;
    p_values = [];
    x = linspace(-1, 1, LF);
    
    % ������������ �������. 
    y = (1./(x.^2+1))';
    % ====================
   
    C = zeros(LF,N);
     for i = 1:LF
         for j = 1:N
             C(i,j) = x(i)^(j-1);
         end
     end
    WC = C;

    a = C\y; 
    A = C*a;
    pk = 2; 
    % k - ����������� �� ���������� ��������. 
    for k = 1:250
        iter = iter + 1;
        pk = min([p, K*pk]); % �������� ����������
        p_values = [p_values; pk];
        if isnan(A - y)
            disp([A]);  % => NAN 
            disp([a]); % => NAN 
            disp([w]);
            disp([w1]); % => � ���� 
            break
        end
        e = A - y; % ������ ������
        
        % ����� ����������� ������
        [ymax,imax,ymin,imin] = extrema(e);
        if max(ymax) - min(ymax) > (max(abs(ymin)) - min(abs(ymin)))
            mmErrors = [mmErrors;  max(ymax) - min(ymax)];
        else
            mmErrors = [mmErrors; (max(abs(ymin)) - min(abs(ymin)))];
        end
        
        % ���������� �����
        w1 = abs(e.^((pk-2)/2)); 
        w = w1/sum(w1);
        
        % �������� �� ���������
        if isnan(w)
            disp(['nan']);
            break
        end

        for m = 1:N
            WC(:,m) = w.*C(:,m); 
        end
        a = WC\(w.*y); 
        q = 1/(pk-1); % �������� �����
        A = q*C*a + (1-q)*A; % ��������� ���������� �����

        Errors = [Errors; max(abs(e))];
        if max(abs(e)) < 1e-15
           vpa(e)
           break
        end
        
    end

    % ������ ������������ � ������������ �������.
    figure(1)
    plot(x, A, x, y), set(gca,'FontSize',18)
    legend('r(x)', 'f(x)'), grid on
    % ������ ������� ������. 
    figure(2)
    plot(x,e, 'Linewidth',2), set(gca,'FontSize',18)
    legend('Error function'),grid on
end

[min_err, bestN] = min(Errors);

% ����� ������������� ���������� � �������� p
min_err
vpa(pk)


    
    


