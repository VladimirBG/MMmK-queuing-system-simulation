function MMmKim
clc, close, clear

%% Задание значений начальных параметров СМО
DialogM = inputdlg({'Средняя интенсивность потока поступления заявок (Lambda)' 'Средняя интенсивность потока обслуживания заявок (Mu)' 'Максимальная длина очереди (m)' 'Допустимое число требований в системе (K)' 'Время моделирования СМО' 'Число повторов имитации'}, 'Введите начальные параметры системы', [1; 1; 1; 1; 1; 1], {'1.254' '0.756' '3' '5' '20' '50'});

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

    %% Цикл для каждого момента времени treq
    for i = 1 : length(treq)

        while tcur < treq(i)
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
fprintf(' Система моделировалась N = %i раз\n', N);
fprintf('\n\t СРЕДНИЕ ЗНАЧЕНИЯ РЕЗУЛЬТАТОВ МОДЕЛИРОВАНИЙ:\n');
fprintf(' Поступило заявок Creq = %f\n', CreqAve/(N-1));
fprintf(' Обслужено заявок Cserv = %f\n', CservAve/N);
fprintf(' Заявок отклонено Creq = %f\n', CrejectedAve/N);
fprintf(' Заявок в очереди на момент окончания моделирования QueueAve = %f\n', QueueAve/N);
fprintf(' Заявок, в процессе обслуживания, на момент окончания моделирования DsAve = %f\n', DsAve/N);

