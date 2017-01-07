function m2
clc, clear

L = 4; % Интенсивность поступления заявок
M = 1; % Интенсивнось обслуживания заявок
m = 2; % Максимальная очередь заявок
K = 5; % Параметр К
Tm = 100; % Время имитационного моделирования СМО
N = 50; % Число повторов моделирования


 
    %% 
CreqAve = 0; % Счетчик среднего числа поступивших заявок
CservAve = 0; % Счетчик среднего числа обслуженных заявок
CrejectedAve = 0; % Счетчик среднего числа отклоненных заявок
QueueAve = 0; % Счетчик среднего числа заявок в очереди на момент окончания моделирования
DsAve = 0; % Счетчик среднего числа заявок, находящиеся в процессе обслуживания, на момент окончания моделирования
PdsBusy = zeros(1, K-m+2); % Вероятность загрузки приборов обслуживания (K-m+1 - Все приборы свободны, K-m+2 - Все приборы заняты)
PdsBusyNotAll = 0; % Вероятность загрузки одного прибора обслуживания
QsCount = zeros(1, m+1); % Средние значения нахождения в очереди i заявок

for I = 1 : N
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
    Dstatus = zeros(1, K-m); % Статус прибора (0-свободен, !0 - занят)
    Queue = 0; % Число заявок в очереди
    treq = 0; % Контейнер с моментами поступления заявок
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
        %Освободились ли
        for J = 1:length(Dstatus)
            if Dstatus(J) > 0 && Dstatus(J) < CT
                Cserv = Cserv + 1;
                Dstatus(J) = 0;
            end
        end
        
        %Обработка очереди
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
                %Заявка епты
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
fprintf(' Вероятность отказа Pnot = %f\n', Pnot);