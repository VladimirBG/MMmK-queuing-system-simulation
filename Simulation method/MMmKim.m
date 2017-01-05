function MMmKim
clc, close, clear

%% ������� �������� ��������� ���������� ���
DialogM = inputdlg({'������� ������������� ������ ����������� ������ (Lambda)' '������� ������������� ������ ������������ ������ (Mu)' '������������ ����� ������� (m)' '���������� ����� ���������� � ������� (K)' '����� ������������� ���' '����� �������� ��������'}, '������� ��������� ��������� �������', [1; 1; 1; 1; 1; 1], {'1.254' '0.756' '3' '5' '20' '50'});

L = str2double(DialogM{1}); % ������������� ����������� ������
M = str2double(DialogM{2}); % ������������ ������������ ������
m = str2double(DialogM{3}); % ������������ ������� ������
K = str2double(DialogM{4}); % �������� �
Tm = str2double(DialogM{5}); % ����� ������������� ������������� ���
N = str2double(DialogM{6}); % ����� �������� �������������

%% 
CreqAve = 0; % ������� �������� ����� ����������� ������
CservAve = 0; % ������� �������� ����� ����������� ������
CrejectedAve = 0; % ������� �������� ����� ����������� ������
QueueAve = 0; % ������� �������� ����� ������ � ������� �� ������ ��������� �������������
DsAve = 0; % ������� �������� ����� ������, ����������� � �������� ������������, �� ������ ��������� �������������

for Nt = 1 : N % ����� �������� �������������
        %% 
    tcur = 0; % ������� ����� ������ �������
    tcurCorrected = 0; % ��������������� ������� �����
    Dstatus = zeros(1, K-m); % ������ ������� (0-��������, !0 - �����)
    Queue = 0; % ����� ������ � �������
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
                tserv = -1/M*log(rand);
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
                Queue = Queue + 1;
            end

            if busyFlag && Queue > m %���� ��������� �������
                Crejected = Crejected + 1;
                Queue = Queue - 1;
            else
                if Queue > 0 && not(busyFlag) % ���� ���� �������
                    for qu = 1 : Queue
                        tserv = -1/M*log(rand);
                        for ds = 1 : length(Dstatus)
                            if Dstatus(ds) == 0
                                Dstatus(ds) = tcur + tserv;
                                Queue = Queue - 1;
                                break;
                            end
                        end 
                    end
                end
            end 
        end
    end
    
    CreqAve = CreqAve + Creq; 
    CservAve = CservAve + Cserv; 
    CrejectedAve = CrejectedAve + Crejected;
    QueueAve = QueueAve + Queue;
    DsAve = DsAve + length(Dstatus);
end
fprintf(' ������� �������������� N = %i ���\n', N);
fprintf('\n\t ������� �������� ����������� �������������:\n');
fprintf(' ��������� ������ Creq = %f\n', CreqAve/(N-1));
fprintf(' ��������� ������ Cserv = %f\n', CservAve/N);
fprintf(' ������ ��������� Creq = %f\n', CrejectedAve/N);
fprintf(' ������ � ������� �� ������ ��������� ������������� QueueAve = %f\n', QueueAve/N);
fprintf(' ������, � �������� ������������, �� ������ ��������� ������������� DsAve = %f\n', DsAve/N);

