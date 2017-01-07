function MMmKim
clc, close, clear

%% Задание значений начальных параметров СМО
DialogM = inputdlg({'Средняя интенсивность потока поступления заявок (Lambda)' 'Средняя интенсивность потока обслуживания заявок (Mu)' 'Максимальная длина очереди (m)' 'Допустимое число требований в системе (K)' 'Время моделирования СМО' 'Число повторов имитации'}, 'Введите начальные параметры системы', [1; 1; 1; 1; 1; 1], {'3.52' '0.678' '4' '7' '100' '50'});

L = str2double(DialogM{1}); % Интенсивность поступления заявок
M = str2double(DialogM{2}); % Интенсивнось обслуживания заявок
m = str2double(DialogM{3}); % Максимальная очередь заявок
K = str2double(DialogM{4}); % Параметр К
Tm = str2double(DialogM{5}); % Время имитационного моделирования СМО
N = str2double(DialogM{6}); % Число повторов моделирования

%% 
CreqAve = 0; % Счетчик среднего числа поступивших заявок
CservAve = 0; % Счетчик среднего числа обслуженных заявок
CrejectedAve = 0; % Счетчик среднего числа отклоненных заявок
QueueAve = 0; % Счетчик среднего числа заявок в очереди на момент окончания моделирования
DsAve = 0; % Счетчик среднего числа заявок, находящиеся в процессе обслуживания, на момент окончания моделирования
PdsBusy = zeros(1, K-m+2); % Вероятность загрузки приборов обслуживания (K-m+1 - Все приборы свободны, K-m+2 - Все приборы заняты)
PdsBusyNotAll = 0; % Вероятность загрузки одного прибора обслуживания
QsCount = zeros(1, m+1); % Средние значения нахождения в очереди i заявок



for Nt = 1 : N % Число повторов моделирования
        %% 
    tcur = 0; % Текущее время работы системы
    tcurCorrected = 0; % Корректированое текущее время
    Dstatus = zeros(1, K-m); % Статус прибора (0-свободен, !0 - занят)
    Queue = 0; % Число заявок в очереди
    treq = 0; % Контейнер с моментами поступления заявок
    tserv = 0; % Время обслуживания заявки

    Creq = 0; % Счетчик поступивших заявок
    Cserv = 0; % Счетчик обслуженных заявок
    Crejected = 0; % Счетчик отклоненных заявок
    
    %%  Генерация матрицы моментов поступления заявок
    ttemp = 0;
    while (treq(end) <= Tm)
        Creq = Creq + 1;
        ttemp = ttemp + (-1/L*log(rand));
        treq(Creq) = ttemp;
    end
    Creq = Creq - 1;
    treq(end) = [];

    %% Цикл для каждого момента времени treq
    for i = 1 : length(treq)

        while tcur < treq(i)
             
            % Счетчики наличия очереди и занятости приборов обслуживания
            flDstatus = true; % флаг все приборы свободны
            flDstatusBAll = true; % Флаг занятости всех приборов 
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
      
            
            tcurCorrected = Dstatus(1);  % Ближайшее время
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

            for ds = 1 : length(Dstatus) % Если заявка обработалась, то статус прибора = 0
                if Dstatus(ds) <= tcur && Dstatus(ds) > 0
                    Cserv = Cserv + 1;
                    Dstatus(ds) = 0;
                end
            end

            % Проверка наличия свободных приборов
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

            if busyFlag  % Если все прибоы заняты, то инкрементируется счетчик очереди
                Queue = Queue + 1;
            end

            if busyFlag && Queue > m %Если превысило очередь
                Crejected = Crejected + 1;
                Queue = Queue - 1;
            else
                if Queue > 0 && not(busyFlag) % если есть очередь
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
    
fprintf(' Система моделировалась N = %i раз\n', N);
fprintf('\n\t СРЕДНИЕ ЗНАЧЕНИЯ РЕЗУЛЬТАТОВ МОДЕЛИРОВАНИЯ СМО M/M/%i/%i:\n', m, K);
fprintf(' Поступило заявок Creq = %f\n', CreqAve);
fprintf(' Обслужено заявок Cserv = %f\n', CservAve);
fprintf(' Заявок отклонено Creq = %f\n', CrejectedAve);
fprintf(' Заявок в очереди на момент окончания моделирования QueueAve = %f\n', QueueAve);
fprintf(' Заявок, в процессе обслуживания, на момент окончания моделирования DsAve = %f\n', DsAve);

Pnot = CrejectedAve / (CrejectedAve + CservAve);
fprintf(' Вероятность отказа Pnot = %f\n', Pnot);
Q = 1 - Pnot;
fprintf(' Относительная пропускная способность Q = %f\n', Q);
Ab = L*Q;
fprintf(' Абсолютная пропускная способность A = %f\n', Ab);
QsC = 0;
for Ji = 2 : length(QsCount)
   QsC = QsC + QsCount(Ji); 
end
Pq = QsC / (QsC + QsCount(1));
fprintf(' Вероятность наличия очереди Pq = %f\n', Pq);
Ps = PdsBusy(K-m+2) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1));
fprintf(' Вероятность загрузки всех каналов обслуживания Ps = %f\n', Ps);
Ns = (Creq)/Tm;
fprintf(' Среднее количество требований в системе Ns = %f\n', Ns);
for Ji = 1 : K-m
  PpB = PdsBusy(Ji) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1)); 
  fprintf(' Загрузка %i прибора: %f\n', Ji, PpB);  
end

%% Расчет доверительных интервалов
fprintf('\n\t Доверительные интервалы:\n');

PnotC = vpa(CInterval95(Pnot,N));
fprintf(' Вероятность отказа Pnot:\n   От %f\n   До %f\n', PnotC);
QC = vpa(CInterval95(Q,N));
fprintf(' Относительная пропускная способность Q:\n   От %f\n   До %f\n', 1-PnotC);
AbC = vpa(CInterval95(Ab,N));
fprintf(' Абсолютная пропускная способность A:\n   От %f\n   До %f\n', QC*L);
PqC = vpa(CInterval95(Pq,N));
fprintf(' Вероятность наличия очереди Pq:\n   От %f\n   До %f\n', PqC);
PsC = vpa(CInterval95(Ps,N));
fprintf(' Вероятность загрузки всех каналов обслуживания Ps:\n   От %f\n   До %f\n', PsC);
NsC = vpa(CInterval95(Ns,N));
fprintf(' Среднее количество требований в системе Ns:\n   От %f\n   До %f\n', NsC);
for Ji = 1 : K-m
  PpB = PdsBusy(Ji) / (PdsBusyNotAll + PdsBusy(K-m+2) + PdsBusy(K-m+1));
  PpBC = vpa(CInterval95(PpB,N));
  for I = 1 : length(PpBC)
      if PpBC(I) > 1
          PpBC(I) = 1;
      end
  end
  fprintf(' Загрузка %i прибора::\n   От %f\n   До %f\n', Ji, PpBC);  
end

function f = CInterval95(xc,Ni)
%% Функция расчета доверительных интервалов с p = 0.95
syms tot;
E = {(xc - tot)^2 - 1.96^2*tot/Ni == 0};
f = solve(E{:});
