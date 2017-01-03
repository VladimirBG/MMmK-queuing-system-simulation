function MMmKim
clc, close, clear

%% ������� �������� ��������� ���������� ���
DialogM = inputdlg({'������� ������������� ������ ����������� ������ (Lambda)' '������� ������������� ������ ������������ ������ (Mu)' '������������ ����� ������� (m)' '���������� ����� ���������� � ������� (K)' '����� ������������� ���' '����� �������� ��������'}, '������� ��������� ��������� �������', [1; 1; 1; 1; 1; 1], {'1.254' '0.756' '3' '5' '20' '50'});
DialogT = inputdlg({'�� T1' '�� T2'}, '�������� ������� ������������ �������� [T1 - T2]', [1 85; 1 85], {'5', '10'});

L = str2double(DialogM{1}); % ������������� ����������� ������
M = str2double(DialogM{2}); % ������������ ������������ ������
m = str2double(DialogM{3}); % ������������ ������� ������
K = str2double(DialogM{4}); % �������� �
Tm = str2double(DialogM{5}); % ����� ������������� ������������� ���
N = str2double(DialogM{5}); % ����� �������� �������������

T1 = str2double(DialogT{1}); % ��������� ������� ������������
T2 = str2double(DialogT{2});

%% 

tcur = 0; % ������� ����� ������ �������
tcurCorrected = 0; % ��������������� ������� �����
Dstatus = zeros(1, K-m); % ������ ������� (0-��������, !0 - �����)
Queues = 0; % ����� ������ � �������
treq = 0; % ��������� � ��������� ����������� ������
tserv = 0; % ����� ������������ ������

Creq = 0; % ������� ����������� ������
Cserv = 0; % ������� ����������� ������
Crejected = 0; % ������� ����������� ������

%%  ��������� ������� �������� ����������� ������
rng('shuffle');
ttemp = 0;
while (treq(end) <= Tm)
    Creq = Creq + 1;
    ttemp = ttemp + (-1/L*log(rand));
    treq(Creq) = ttemp;
end
Creq = Creq - 1;
treq(end) = [];
stairs(treq);

%% ���� ��� ������� ������� ������� treq
for i = 1 : length(treq)
    
    while tcur < treq(i)
        tcurCorrected = Dstatus(1);  % ��������� �����
        for ds = 1 : length(Dstatus)
             if not(Dstatus(ds) == 0) && Dstatus(ds) < tcurCorrected
                tcurCorrected = Dstatus(ds);
             end
        end 

        if tcurCorrected < treq(i) && not(tcurCorrected == 0)
            tcur = tcurCorrected;
            flag = false;
        else
            tcur = treq(i);
            flag = true;
        end
          
        for ds = 1 : length(Dstatus) % ���� ������ ������������, �� ������ ������� = 0
            if Dstatus(ds) <= tcur && Dstatus(ds) > 0
                Cserv = Cserv + 1;
                Dstatus(ds) = 0;
            end
        end

        % �������� ������� ��������� ��������
        if flag
            busyFlag = true;
            rng('shuffle');
            tserv = 1/M*log(randi([T1,T2]));
            for ds = 1 : length(Dstatus)
                if Dstatus(ds) == 0
                    Dstatus(ds) = tcur + tserv;
                    busyFlag = false; 
                    break;
                end
            end
        else
            busyFlag = false;
        end

        if busyFlag  % ���� ��� ������ ������, �� ���������������� ������� �������
            Queues = Queues + 1;
        end

        if busyFlag && Queues > m %���� ��������� �������
            Crejected = Crejected + 1;
            Queues = Queues - 1;
        else
            if Queues > 0 && not(busyFlag) % ���� ���� �������
                for qu = 1 : Queues
                    rng('shuffle');
                    tserv = 1/M*log(randi([T1,T2]));
                    for ds = 1 : length(Dstatus)
                        if Dstatus(ds) == 0
                            Dstatus(ds) = tcur + tserv;
                            Queues = Queues - 1;
                            break;
                        end
                    end 
                end
            end
        end 
    end
end

display(Cserv);
display(Crejected);
display(Creq);
display(Dstatus);
display(Queues);
