function mmmk
clc,close 
% ��������� ��������� �������
DialogM = inputdlg({'������� ������������� ������ ����������� ������ (Lambda)' '������� ������������� ������ ������������ ������ (Mu)' '������������ ����� ������� (m)' '���������� ����� ���������� � ������� (K)'}, '������� ��������� ��������� �������', [1; 1; 1; 1], {'1.254' '0.756' '2' '4'});

L = str2double(DialogM{1});
M = str2double(DialogM{2});
m = str2double(DialogM{3});
K = str2double(DialogM{4}); 
%%%%%%%%%%%%%%%%%%%%%%%%
global Pij
Pij = zeros(K+1,K+1)';
Pij(1, 1:2) = [-L, M]; 
I = 1;
while I < K+1
    I = I+1;
    if I == K+1
         Pij(I, I-1) = L;
         Pij(I, I) = -M*(K-m+1);
    elseif I <= K-m+1
        Pij(I, I-1) = L;
        Pij(I, I) = -(L+M*(I-1));
        Pij(I, I+1) = M*I;
        else
            Pij(I, I-1) = L;
            Pij(I, I) = -(L+M*(K-m+1));
            Pij(I, I+1) = M*(K-m+1);
    end   
end

%% ��������� �������������� ����. ���������
P0 = [1;zeros(length(Pij)-1,1)];
T = [0,20];
[t,P] = ode23(@cmo, T, P0);
 
%% ���������� ��������� ������������ ���������
%% plot %% � ���������� �������
plot(t,P, 'LineWidth', 2);

grid on
N = length(Pij)-1;
arr = [0:N]';
str = num2str(arr);
legend(strcat('\bf\itP\rm\bf_', str, '(\itt\rm\bf)'));
title(sprintf('%s ����������� ��������� ������� M/M/%d/%d', '\bf\fontsize{12}',m, K));
xlabel('\bf\it\fontsize{12}  t')
ylabel('\bf\fontsize{12}\itP\rm\bf(\itt\rm\bf)');
set(gca, 'fontweight','bold', 'fontsize',10)

fprintf('\n ������������ �����������:\n');
for J = 1 : length(Pij)
    fprintf('\tP%d = %f\n', J-1, P(end,J));
end

%% ��������� ����������� �������
fprintf('\n ��������� �����������:\n');
A = Pij;
A(K+2,:)=ones(1,length(Pij));
fa = zeros(1,length(Pij))';
fa(K+2,:) = 1;
Pf = vpa(mldivide(A, fa), 7);
for J = 1 : length(Pf)
    fprintf('\tPf%d = %f\n', J-1, Pf(J));
end


fprintf('\n\t ������������ ��������������:\n');
Pnot = P(end,end);
fprintf(' ����������� ������ Pnot = %f\n', P(end,end));
Q = 1 - Pnot;
fprintf(' ������������� ���������� ����������� Q = %f\n', Q);
Ab = L*Q;
fprintf(' ���������� ���������� ����������� A = %f\n', Ab);
Pq = sum(P(end, m+1:end));
fprintf(' ����������� ������� ������� Pq = %f\n', Pq);
Ps = sum(P(end, m:end));
fprintf(' ����������� �������� ���� ������� ������������ Ps = %f\n', Ps);
Ns = [0:length(Pij)-1]*P(end,:)';
fprintf(' ������� ���������� ���������� � ������� Ns = %f\n', Ns);
fprintf(' ������� ����� ���������� ���������� � ������� Ts = %f\n', Ns/L);

function f = cmo(t,P);
%% ������� �������� ������ ������ 
%% ���������������� ���������
global Pij
f = Pij*P;
