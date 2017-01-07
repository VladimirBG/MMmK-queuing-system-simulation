function MMmKim
clc, close, clear

%% ������� �������� ��������� ���������� ���
DialogM = inputdlg({'������� ������������� ������ ����������� ������ (Lambda)' '������� ������������� ������ ������������ ������ (Mu)' '������������ ����� ������� (m)' '���������� ����� ���������� � ������� (K)' '����� ������������� ���' '����� �������� ��������'}, '������� ��������� ��������� �������', [1; 1; 1; 1; 1; 1], {'3.52' '0.678' '4' '7' '100' '50'});

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
PdsBusy = zeros(1, K-m+2); % ����������� �������� �������� ������������ (K-m+1 - ��� ������� ��������, K-m+2 - ��� ������� ������)
PdsBusyNotAll = 0; % ����������� �������� ������ ������� ������������
QsCount = zeros(1, m+1); % ������� �������� ���������� � ������� i ������



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
    ttemp = 0;
    while (treq(end) <= Tm)
        Creq = Creq + 1;
        ttemp = ttemp + (-1/L*log(rand));
        treq(Creq) = ttemp;
    end
    Creq = Creq - 1;
    treq(end) = [];

    %% ���� ��� ������� ������� ������� treq
    for i = 1 : length(treq)

        while tcur < treq(i)
             
            % �������� ������� ������� � ��������� �������� ������������
            flDstatus = true; % ���� ��� ������� ��������
            flDstatusBAll = true; % ���� ��������� ���� �������� 
            for Ji = 1 : length(Dstatus)
                if Dstatus(Ji) > 0
                    PdsBusy(Ji) = PdsBusy(Ji) + 1;
                    flDstatus = flDstatus && false;
                    flDstatusBAll = flDstatusBAll && true;
                else
                    flDstatusBAll = flDstatusBAll && false;
                end
            end
            if flDstatus
                PdsBusy(K - m + 1) = PdsBusy(K - m + 1) + 1;
            elseif flDstatusBAll
                PdsBusy(K - m + 2) = PdsBusy(K - m + 2) + 1;
            else
                PdsBusyNotAll = PdsBusyNotAll + 1;
            end
            QsCount(Queue + 1) = QsCount(Queue + 1) + 1;
      
            
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

CreqAve = CreqAve/N; 
CservAve = CservAve/N; 
CrejectedAve = CrejectedAve/N;
QueueAve = QueueAve/N;
DsAve = DsAve/N;
QsCount = QsCount ./ N;
PdsBusy = PdsBusy ./ N;
PdsBusyNotAll = PdsBusyNotAll / N;
    
fprintf(' ������� �������������� N = %i ���\n', N);
fprintf('\n\t ������� �������� ����������� ������������� ��� M/M/%i/%i:\n', m, K);
fprintf(' ��������� ������ Creq = %f\n', CreqAve);
fprintf(' ��������� ������ Cserv = %f\n', CservAve);
fprintf(' ������ ��������� Creq = %f\n', CrejectedAve);
fprintf(' ������ � ������� �� ������ ��������� ������������� QueueAve = %f\n', QueueAve);
fprintf(' ������, � �������� ������������, �� ������ ��������� ������������� DsAve = %f\n', DsAve);

Pnot = CrejectedAve / (CrejectedAve + CservAve);
fprintf(' ����������� ������ Pnot = %f\n', Pnot);
Q = 1 - Pnot;
fprintf(' ������������� ���������� ����������� Q = %f\n', Q);
Ab = L*Q;
fprintf(' ���������� ���������� ����������� A = %f\n', Ab);
QsC = 0;
for Ji = 2 : length(QsCount)
   QsC = QsC + QsCount(Ji); 
end
Pq = QsC / (QsC + QsCount(1));
fprintf(' ����������� ������� ������� Pq = %f\n', Pq);
Ps = PdsBusy(K-m+2) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1));
fprintf(' ����������� �������� ���� ������� ������������ Ps = %f\n', Ps);
Ns = (Creq)/Tm;
fprintf(' ������� ���������� ���������� � ������� Ns = %f\n', Ns);
for Ji = 1 : K-m
  PpB = PdsBusy(Ji) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1)); 
  fprintf(' �������� %i �������: %f\n', Ji, PpB);  
end

%% ������ ������������� ����������
fprintf('\n\t ������������� ���������:\n');

PnotC = vpa(CInterval95(Pnot,N));
fprintf(' ����������� ������ Pnot:\n   �� %f\n   �� %f\n', PnotC);
QC = vpa(CInterval95(Q,N));
fprintf(' ������������� ���������� ����������� Q:\n   �� %f\n   �� %f\n', 1-PnotC);
AbC = vpa(CInterval95(Ab,N));
fprintf(' ���������� ���������� ����������� A:\n   �� %f\n   �� %f\n', QC*L);
PqC = vpa(CInterval95(Pq,N));
fprintf(' ����������� ������� ������� Pq:\n   �� %f\n   �� %f\n', PqC);
PsC = vpa(CInterval95(Ps,N));
fprintf(' ����������� �������� ���� ������� ������������ Ps:\n   �� %f\n   �� %f\n', PsC);
NsC = vpa(CInterval95(Ns,N));
fprintf(' ������� ���������� ���������� � ������� Ns:\n   �� %f\n   �� %f\n', NsC);
for Ji = 1 : K-m
  PpB = PdsBusy(Ji) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1));
  PpBC = vpa(CInterval95(PpB,N));
  for I = 1 : length(PpBC)
      if PpBC(I) > 1
          PpBC(I) = 1;
      end
  end
  fprintf(' �������� %i �������::\n   �� %f\n   �� %f\n', Ji, PpBC);  
end

function f = CInterval95(xc,Ni)
%% ������� ������� ������������� ���������� � p = 0.95
syms tot;
E = {(xc - tot)^2 - 1.96^2*tot/Ni == 0};
f = solve(E{:});
