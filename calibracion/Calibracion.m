% =========================================================================
% CALIBRACIÓN DE LOS SENSORES
% =========================================================================
% En este archivo de MATLAB se proporciona el código para poder calibrar
% tanto los sensores propioceptivos como los sensores exteroceptivos del
% robot. Para ello carga el entorno llamado "calibracion.xml" y ejecuta
% el código.

%% Calibración de los sensores propioceptivos - Desplazamiento

clear all;
close all;
clc;

v_x = 0.2; % Velocidad lineal en m/s
t = 10; % Tiempo total de avance
t_res = 0.1; % Time-step del avance
matrix_errors = []; % Matriz para guardar en cada fila los errores

for j = 1:100
    % Reseteo de las condiciones iniciales para comenzar la calibración
    apoloPlaceMRobot('Camarero', [0 0 0], 0);
    apoloResetOdometry('Camarero');
    apoloUpdate();
    
    % Calibración odometría desplazamiento en línea recta
    for i = 1:(t/t_res)
        apoloMoveMRobot('Camarero', [v_x, 0], t_res);
        apoloUpdate()
        odometry = apoloGetOdometry('Camarero');
        real_location = apoloGetLocationMRobot('Camarero');
        matrix_errors = [matrix_errors ; odometry(1) - real_location(1)];
        apoloResetOdometry ('Camarero', [real_location(1), real_location(2), real_location(4)]);
    end
end

% Calcular media y varianza
mean_val = mean(matrix_errors);
variance_val = var(matrix_errors);

% Ajuste a una distribución normal
x = linspace(min(matrix_errors), max(matrix_errors), 100);
normal_curve = (1/sqrt(2*pi*variance_val))*exp(-0.5*((x - mean_val).^2 / variance_val));

% Dibujar el histograma y la curva normal de ajuste
figure;
histogram(matrix_errors, 24);
hold on;
plot(x, normal_curve * (numel(matrix_errors) * (max(matrix_errors) - min(matrix_errors)) / 24), 'r', 'LineWidth', 2);
xlabel('Error en el desplazamiento');
ylabel('Frecuencia');
title('Histograma de Calibración de Posición');
saveas(gcf, 'calibration_position.epsc', 'epsc');

%% Calibración de los sensores propioceptivos - Giro

clear all;
close all;
clc;

v_x = 0.2; % Velocidad lineal en m/s
R = 2; % Radio de la circunferencia en m
omega = v_x/R; % Velocidad angular en rad/s  
n = 10; % Número de vueltas
t_res = 0.1; % Time-step del giro
n_iteraciones = ceil(n * 2 * pi/(omega * t_res));
matrix_errors = []; % Matriz para guardar en cada fila los errores

% Reseteo de las condiciones iniciales para comenzar la calibración
apoloPlaceMRobot('Camarero', [0 -R 0], 0);

for i = 1:n_iteraciones
    apoloMoveMRobot('Camarero', [v_x, omega], t_res);
    apoloUpdate()
    odometry = apoloGetOdometry('Camarero');
    real_location = apoloGetLocationMRobot('Camarero');
    matrix_errors = [matrix_errors ; odometry(3) - real_location(4)];
    apoloResetOdometry ('Camarero', [real_location(1), real_location(2), real_location(4)]);
end

% Calcular media y varianza
mean_val = mean(matrix_errors);
variance_val = var(matrix_errors);

% Ajuste a una distribución normal
x = linspace(min(matrix_errors), max(matrix_errors), 100);
normal_curve = (1/sqrt(2*pi*variance_val))*exp(-0.5*((x - mean_val).^2 / variance_val));

% Dibujar el histograma y la curva normal de ajuste
figure;
histogram(matrix_errors, 24);
hold on;
plot(x, normal_curve * (numel(matrix_errors) * (max(matrix_errors) - min(matrix_errors)) / 24), 'r', 'LineWidth', 2);
xlabel('Error en el giro');
ylabel('Frecuencia');
title('Histograma de Calibración del Giro');
saveas(gcf, 'calibration_rotation.epsc', 'epsc');

%% Calibración de los sensores exteroceptivos - Ultrasonidos

clear all;
close all;
clc;

n_iteraciones = 10000;
matrix_errors = []; % Matriz para guardar en cada fila los errores

% Reseteo de las condiciones iniciales para comenzar la calibración
apoloPlaceMRobot('Camarero', [3 0 0], 0);
apoloResetOdometry('Camarero');
apoloUpdate();

for i = 1:n_iteraciones
    ultrasonic = apoloGetUltrasonicSensor('uc0');
    real_location = apoloGetLocationMRobot('Camarero');
    matrix_errors = [matrix_errors ; ultrasonic - (3.8 - real_location(1))];
    apoloResetOdometry ('Camarero', [real_location(1), real_location(2), real_location(4)]);
end

% Calcular media y varianza
mean_val = mean(matrix_errors);
variance_val = var(matrix_errors);

% Ajuste a una distribución normal
x = linspace(min(matrix_errors), max(matrix_errors), 100);
normal_curve = (1/sqrt(2*pi*variance_val))*exp(-0.5*((x - mean_val).^2 / variance_val));

% Dibujar el histograma y la curva normal de ajuste
figure;
histogram(matrix_errors, 100);
hold on;
plot(x, normal_curve * (numel(matrix_errors) * (max(matrix_errors) - min(matrix_errors)) / 100), 'r', 'LineWidth', 2);
xlabel('Error en la medida');
ylabel('Frecuencia');
title('Histograma de Calibración de los Sensores de Ultrasonidos');
saveas(gcf, 'calibration_ultrasonic.epsc', 'epsc');

%% Calibración de los sensores exteroceptivos - Telémetro Láser LMS100

clear all;
close all;
clc;

n_iteraciones = 10000; % Número de medidas del láser de la baliza
matrix_errors_a = []; % Matriz para guardar en cada fila los errores ángulo
matrix_errors_d = []; % Matriz para guardar en cada fila los errores distancia

% Reseteo de las condiciones iniciales para comenzar la calibración
apoloPlaceMRobot('Camarero', [0 0 0], 0);
apoloResetOdometry('Camarero');
apoloUpdate();

% Tomar medidas n veces
for i = 1:n_iteraciones
    laser = apoloGetLaserLandMarks('LMS100');
    matrix_errors_a = [matrix_errors_a ; laser.angle - atan(2/3)];
    matrix_errors_d = [matrix_errors_d ; laser.distance - sqrt(13)];
end

% Calcular medias y varianzas
mean_val_a = mean(matrix_errors_a);
mean_val_d = mean(matrix_errors_d);
variance_val_a = var(matrix_errors_a);
variance_val_d = var(matrix_errors_d);

% Ajustes a distribuciones normales
x_a = linspace(min(matrix_errors_a), max(matrix_errors_a), 100);
x_d = linspace(min(matrix_errors_d), max(matrix_errors_d), 100);
normal_curve_a = (1/sqrt(2*pi*variance_val_a))*exp(-0.5*((x_a - mean_val_a).^2 / variance_val_a));
normal_curve_d = (1/sqrt(2*pi*variance_val_d))*exp(-0.5*((x_d - mean_val_d).^2 / variance_val_d));

% Dibujar los histogramas y las curvas normales de ajuste
figure;
subplot(1,2,1);
    histogram(matrix_errors_a, 24);
    hold on;
    plot(x_a, normal_curve_a * (numel(matrix_errors_a) * (max(matrix_errors_a) - min(matrix_errors_a)) / 24), 'r', 'LineWidth', 2);
    xlabel('Error en la medida');
    ylabel('Frecuencia');
    title('Histograma de Calibración del Telémetro Láser - Ángulo');
subplot(1,2,2);
    histogram(matrix_errors_d, 24);
    hold on;
    plot(x_d, normal_curve_d * (numel(matrix_errors_d) * (max(matrix_errors_d) - min(matrix_errors_d)) / 24), 'r', 'LineWidth', 2);
    xlabel('Error en la medida');
    ylabel('Frecuencia');
    title('Histograma de Calibración del Telémetro Láser - Distancia');
    saveas(gcf, 'calibration_laser.epsc', 'epsc');