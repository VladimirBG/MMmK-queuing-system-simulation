function MMmKim
clc, close, clear

%% Задание значений начальных параметров СМО
DialogM = inputdlg({'Средняя интенсивность потока поступления заявок (Lambda)' 'Средняя интенсивность потока обслуживания заявок (Mu)' 'Максимальная длина очереди (m)' 'Допустимое число требований в системе (K)' 'Время моделирования СМО' 'Число повторов имитации'}, 'Введите начальные параметры системы', [1; 1; 1; 1; 1; 1], {'1.254' '0.756' '3' '5' '60' '50'});
DialogT = inputdlg({'От T1' 'До T2'}, 'Интервал времени обслуживания приборов [T1 - T2]', [1 85; 1 85], {'5', '10'});

L = str2num(DialogM{1}); % Интенсивность поступления заявок
M = str2num(DialogM{2}); % Интенсивнось обслуживания заявок
m = str2num(DialogM{3}); % Максимальная очередь заявок
K = str2num(DialogM{4}); % Параметр К
Tm = str2num(DialogM{5}); % Время имитационного моделирования СМО
N = str2num(DialogM{5}); % Число повторов моделирования

T1 = str2num(DialogT{1}); % Интервалы времени обслуживания
T2 = str2num(DialogT{2});

%% 

tcur = 0; % Текущее время работы системы
tcurCorrected = inf; % Корректированое текущее время
Dstatus = zeros(1, K-m); % Статус прибора (0-свободен, !0 - занят)
Queues = 0; % Число заявок в очереди
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
treq(end) = [];
stairs(treq);

%% Цикл для каждого момента времени treq
for i = 1 : length(treq)
    tcur = treq(i);
    for ds = 1 : length(Dstatus) % Если заявка обработалась, то статус прибора = 0
        if not(Dstatus(ds) > tcur)
            if not(Dstatus(ds) == 0)
                Cserv = Cserv + 1;
            end
            Dstatus(ds) = 0;
        end
    end
    
    % Проверка наличия свободных приборов
    busyFlag = true;
    rng('shuffle');
    tserv = randi([T1,T2]);
    for ds = 1 : length(Dstatus)
        if Dstatus(ds) == 0
            Dstatus(ds) = tcur + tserv;
            busyFlag = false; 
            break;
        end
    end
    
    if busyFlag  % Если все прибоы заняты, то инкрементируется счетчик очереди
        Queues = Queues + 1;
    end
    
    if busyFlag && Queues > m %Если превысило очередь
        Crejected = Crejected + 1;
        Queues = Queues - 1;
    else
        if Queues > 0  % если есть очередь
            for qu = 1 : Queues
                rng('shuffle');
                tserv = randi([T1,T2]);
                for ds = 1 : length(Dstatus)
                    if Dstatus(ds) == 0
                        Dstatus(ds) = tcur + tserv;
                        break;
                    end
                end 
            end
        end
    end
    
    tcurCorrected = Dstatus(1);  % Ближайшее время
    for ds = 1 : length(Dstatus)
         if not(Dstatus(ds) == 0)
            if Dstatus(ds) < tcurCorrected
                tcurCorrected = Dstatus(ds);
            end
         end
     end 
    
    
    
    
    display(Dstatus);
end
