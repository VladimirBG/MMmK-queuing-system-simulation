function m2
clc, clear

L = 4; % ������������� ����������� ������
M = 1; % ������������ ������������ ������
m = 2; % ������������ ������� ������
K = 5; % �������� �
Tm = 100; % ����� ������������� ������������� ���
N = 50; % ����� �������� �������������


 
    %% 
CreqAve = 0; % ������� �������� ����� ����������� ������
CservAve = 0; % ������� �������� ����� ����������� ������
CrejectedAve = 0; % ������� �������� ����� ����������� ������
QueueAve = 0; % ������� �������� ����� ������ � ������� �� ������ ��������� �������������
DsAve = 0; % ������� �������� ����� ������, ����������� � �������� ������������, �� ������ ��������� �������������
PdsBusy = zeros(1, K-m+2); % ����������� �������� �������� ������������ (K-m+1 - ��� ������� ��������, K-m+2 - ��� ������� ������)
PdsBusyNotAll = 0; % ����������� �������� ������ ������� ������������
QsCount = zeros(1, m+1); % ������� �������� ���������� � ������� i ������

for I = 1 : N
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
    ttemp = 0;
    Dstatus = zeros(1, K-m); % ������ ������� (0-��������, !0 - �����)
    Queue = 0; % ����� ������ � �������
    treq = 0; % ��������� � ��������� ����������� ������
    while (treq(end) <= Tm)
        Creq = Creq + 1;
        ttemp = ttemp + (-1/L*log(rand));
        treq(Creq) = ttemp;
    end
    Creq = Creq - 1;
    treq(end) = [];
      treq = sort(treq, 'descend');
    
    tqq = 0.0001;
    
    for CT = 0:tqq:Tm
        %������������ ��
        for J = 1:length(Dstatus)
            if Dstatus(J) > 0 && Dstatus(J) < CT
                Cserv = Cserv + 1;
                Dstatus(J) = 0;
            end
        end
        
        %��������� �������
        if Queue > 0
            for J = 1:length(Dstatus)
                if Dstatus(J) == 0
                    Queue = Queue - 1;
                    tserv = -1/M*log(rand);
                    Dstatus(J) = tserv + CT;
                end
            end
        end
        
        if ~isempty(treq)
            if(treq(end) <= CT)
                treq(end) = [];
                %������ ����
                tserv = -1/M*log(rand);
                bflag = true;
                for J = 1 : length(Dstatus)
                    if Dstatus(J) == 0
                        Dstatus(J) = tserv + CT;
                        bflag = false;
                        break;
                    end
                end
                if bflag
                    if Queue == m
                        Crejected = Crejected + 1;
                    else
                        Queue = Queue + 1;
                    end
                end  
            end
        end
    end
    CreqAve = CreqAve + Creq;
    CservAve = CservAve + Cserv;
    CrejectedAve = CrejectedAve + Crejected;
end
    CreqAve = CreqAve /N;
    CservAve = CservAve /N;
    CrejectedAve = CrejectedAve/N;
display(CreqAve);
display(CrejectedAve);
display(CservAve);
display(length(Dstatus));
display(Queue);
Pnot = CrejectedAve / CreqAve;
fprintf(' ����������� ������ Pnot = %f\n', Pnot);