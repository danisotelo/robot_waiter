% Inicialización de la posición inicial y la covarianza
Xk = [5; 6; pi/2]; % Posición inicial del robot [x, y, theta]
apoloPlaceMRobot('Camarero', [5 6 0], pi/2);
apoloResetOdometry('Camarero');
apoloUpdate();
Pk = eye(3)*0.000001; % Matriz de covarianza inicial

% Definición de las posiciones de las balizas del restaurante
balizas = [8.7, 0; 8.7, 5; 2.6, 5.5; 23, 5.5; 23, 10; 4, 10; 2.6, 15; 2.6, 17; 21, 17; 21, 15; 21, 10.5; 14, 10];

% Inicialización de las matrices de covarianza del proceso (Qk) y de medida (Rk)
Qk_1 = diag([4.7190e-7 2.4654e-8]); % Suponiendo un pequeño error en el movimiento
R1 = 1.7835e-4; % Varianza en la medida de distancia de una baliza
R2 = 3.6635e-4; % Varianza en la medida del ángulo en radianes


% Variables para almacenar la trayectoria real y estimada
trayectoriaReal = [];
trayectoriaEstimada = [];
Pk_list = {};  % Crear una celda para almacenar Pk en cada iteración
numObservacionesImprobables = 0; % Para contar cuantas observaciones improbables han habido según la distancia de Mahalanobis

% Bucle principal del algoritmo
numIteraciones = 1000; % Define cuántas iteraciones quieres realizar
for k = 1:numIteraciones
    % Movimiento aleatorio del robot
    velLinear = rand * 0.2; % Velocidad lineal aleatoria 
    %velLinear = -rand * 0.2; % Velocidad lineal aleatoria (marcha atrás)
    velAngular = (rand - 0.5) * (pi/10); % Velocidad angular aleatoria 
    timestep = 0.1; % Intervalo de tiempo entre movimientos (100 ms)

    apoloMoveMRobot('Camarero', [velLinear, velAngular], timestep);
    %disp('moving robot')
    apoloUpdate();

    % Obtención de la posición real del robot con apoloGetLocation
    pos = apoloGetLocation('Camarero'); % Esto devuelve un vector de 6 elementos
    
    % Seleccionamos solo los valores de interés: x, y y yaw
    posReal = [pos(1); pos(2); pos(6)]; % Ahora posReal tendrá solo x, y, yaw

    % Obtención de las observaciones de las balizas visibles
    laserLandmarks = apoloGetLaserLandMarks('LMS100');

    % Obtención de Uk con apoloGetOdometry
    odometry = apoloGetOdometry('Camarero');
    disp(odometry)

    % Verificar si odometry(1) 
    if odometry(1) < 0
        % Si odometry(1) negativo, calcular Uk con signo negativo en el desplazamiento
        Uk = [sqrt(odometry(1)^2 + odometry(2)^2); -(2*pi-odometry(3))];
    else
        % Si no, calcular Uk normalmente
        Uk = [sqrt(odometry(1)^2 + odometry(2)^2); odometry(3)];
    end

    if isempty(laserLandmarks.id)
        % No hay balizas visibles, utilizar solo odometría
        disp('No hay balizas visibles, se utiliza la odometría para la predicción del estado.');
       
        % Nuevo ciclo, k-1 = k
        % Actualización de Xk y Pk para el nuevo ciclo
        Xk_1 = Xk;
        Pk_1 = Pk;

        % Predicción del estado
        X_k = [Xk_1(1) + Uk(1)*cos(Xk_1(3) + (Uk(2)/2));
               Xk_1(2) + Uk(1)*sin(Xk_1(3) + (Uk(2)/2));
               Xk_1(3) + Uk(2)];

        % Matriz de transición de estado (Jacobiana del modelo de movimiento
        % respecto al estado)
        Ak = [1 0 (-Uk(1)*sin(Xk_1(3) + Uk(2)/2));
              0 1 (Uk(1)*cos(Xk_1(3) + Uk(2)/2));
              0 0 1];
        % Matriz de control (Jacobiana del modelo de movimiento respecto a la
        % acción de control)
        Bk = [cos(Xk_1(3) + Uk(2)/2) (-0.5*Uk(1)*sin(Xk_1(3) + Uk(2)/2));
              sin(Xk_1(3) + Uk(2)/2) (0.5*Uk(1)*cos(Xk_1(3) + Uk(2)/2));
              0 1];
        
        % Predicción de la covarianza (en base a información de k-1)
        P_k = Ak * Pk_1 * ((Ak)') + Bk * Qk_1 * ((Bk)');
        
        % Solo odometría:
        Xk = X_k;
        Pk = P_k;

        % Almacenar la trayectoria real y la estimada para visualización posterior
        trayectoriaReal = [trayectoriaReal; posReal'];
        trayectoriaEstimada = [trayectoriaEstimada; Xk'];

        % Guardar Pk en la lista para visualización
        Pk_list{k} = Pk;

        apoloResetOdometry('Camarero');
    else
        Zk = []; % Inicialización del vector de observaciones
        visibleIds = []; % Identificadores de balizas visibles
    
        for i = 1:length(laserLandmarks.id)
            Zk = [Zk; laserLandmarks.distance(i); laserLandmarks.angle(i)];
            visibleIds = [visibleIds; laserLandmarks.id(i)];
        end
    
        
        % Nuevo ciclo, k-1 = k
        % Actualización de Xk y Pk para el nuevo ciclo
        Xk_1 = Xk;
        Pk_1 = Pk;
    
        % Predicción del estado
        X_k = [Xk_1(1) + Uk(1)*cos(Xk_1(3) + (Uk(2)/2));
               Xk_1(2) + Uk(1)*sin(Xk_1(3) + (Uk(2)/2));
               Xk_1(3) + Uk(2)];
        
        % Matriz de transición de estado (Jacobiana del modelo de movimiento
        % respecto al estado)
        Ak = [1 0 (-Uk(1)*sin(Xk_1(3) + Uk(2)/2));
              0 1 (Uk(1)*cos(Xk_1(3) + Uk(2)/2));
              0 0 1];
        % Matriz de control (Jacobiana del modelo de movimiento respecto a la
        % acción de control)
        Bk = [cos(Xk_1(3) + Uk(2)/2) (-0.5*Uk(1)*sin(Xk_1(3) + Uk(2)/2));
              sin(Xk_1(3) + Uk(2)/2) (0.5*Uk(1)*cos(Xk_1(3) + Uk(2)/2));
              0 1];
        
        % Predicción de la covarianza (en base a información de k-1)
        P_k = Ak * Pk_1 * ((Ak)') + Bk * Qk_1 * ((Bk)');
    
        % Predicción de la medida Zk_ y cálculo de Hk (usando X_k predicho)
        Zk_ = []; % Inicialización del vector de predicción de observaciones
        Hk = []; % Inicialización de la matriz Jacobiana Hk
    
        for i = 1:length(visibleIds)
            baliza_id = visibleIds(i);
            x_baliza = balizas(baliza_id, 1);
            y_baliza = balizas(baliza_id, 2);
    
            % Distancia y ángulo al i-ésimo landmark visible (función h
            % adaptada al número de balizas visibles)
            distancia = sqrt((X_k(1) - x_baliza)^2 + (X_k(2) - y_baliza)^2);
            angulo = atan2(y_baliza - X_k(2), x_baliza - X_k(1)) - X_k(3);
    
            Zk_ = [Zk_; distancia; angulo]; % Añadir la predicción al vector de observaciones
    
            % Jacobiana para la baliza i (matriz Hk adaptada al número de
            % balizas visibles)
            idx_distancia = 2*i-1; % Índice para la distancia
            idx_angulo = 2*i;      % Índice para el ángulo
            % Derivada parcial respecto X_k(1), X_k(2) y X_k(3) de h
            Hk(idx_distancia:idx_angulo, :) = ...
                [(X_k(1) - x_baliza)/distancia, (X_k(2) - y_baliza)/distancia, 0;
                 (y_baliza - X_k(2))/distancia^2, -(x_baliza - X_k(1))/distancia^2, -1];
        end
    
        % Actualizar Rk basado en las balizas visibles
        Rk_blocks = repmat({diag([R1, R2])}, 1, length(visibleIds));
        Rk = blkdiag(Rk_blocks{:});
        
        % Actualización del estado y la covarianza
        Yk = Zk - Zk_; % Residuo entre la medición real y la predicción
        % Normalizar los ángulos en Yk al rango [-pi, pi]
        for i = 1:length(Yk)/2
            idx_angulo = 2*i;
            if Yk(idx_angulo) > pi
                Yk(idx_angulo) = Yk(idx_angulo) - 2*pi;
            elseif Yk(idx_angulo) < -pi
                Yk(idx_angulo) = Yk(idx_angulo) + 2*pi;
            end
        end
        Sk = Hk * P_k * Hk' + Rk; % Covarianza residual
        Wk = P_k * Hk' * inv(Sk); % Ganancia del filtro de Kalman
        %Xk = X_k + Wk * Yk; % Actualización del estado
        %Pk = (eye(3) - Wk * Hk) * P_k; % Actualización de la covarianza
        
        % Validación del filtro usando la distancia de Mahalanobis
        mahalanobisDistances = zeros(length(visibleIds), 1);
        for i = 1:length(visibleIds)
            idx_distancia = 2*i-1;
            idx_angulo = 2*i;
            Yk_i = Yk(idx_distancia:idx_angulo);
            Sk_i = Sk(idx_distancia:idx_angulo, idx_distancia:idx_angulo);
            mahalanobisDistances(i) = Yk_i' * inv(Sk_i) * Yk_i;
        end
        
        chi2_critical = chi2inv(0.95, 3); % 95% confianza para 3 grados de libertad (x, y, theta)
        disp(mahalanobisDistances)
        improbables = mahalanobisDistances > chi2_critical;
        numObservacionesImprobables = numObservacionesImprobables + sum(improbables);
        
        %if any(improbables)
            %disp(['Observaciones improbables en esta iteración: ', num2str(sum(improbables))]);
        %end
        
        % Identificar las observaciones probables
        observacionesProbables = mahalanobisDistances <= chi2_critical;

        % Filtrar las observaciones improbables
        Yk_probable = Yk(observacionesProbables);
        Hk_probable = Hk(observacionesProbables, :);
        Rk_probable = Rk(observacionesProbables, observacionesProbables);
        
        % Si hay observaciones probables, actualizar el estado y la covarianza
        if ~isempty(Yk_probable)
            Sk = Hk_probable * P_k * Hk_probable' + Rk_probable;
            Wk = P_k * Hk_probable' * inv(Sk);
            Xk = X_k + Wk * Yk_probable;
            Pk = (eye(3) - Wk * Hk_probable) * P_k;
        else
            % Si todas las observaciones son improbables, no actualizamos el estado
            disp('Todas las observaciones son improbables')
            Xk = X_k;
            Pk = P_k;
        end

        % Almacenar la trayectoria real y la estimada para visualización posterior
        trayectoriaReal = [trayectoriaReal; posReal'];
        trayectoriaEstimada = [trayectoriaEstimada; Xk'];

        % Guardar Pk en la lista para visualización
        Pk_list{k} = Pk;

        apoloResetOdometry('Camarero');
    end
end

% Número total de observaciones improbables
disp(['Número total de observaciones improbables: ', num2str(numObservacionesImprobables)]);

% Calcular el error absoluto de la trayectoria (ATE)
diferenciasPosiciones = trayectoriaEstimada(:, 1:2) - trayectoriaReal(:, 1:2);
erroresCuadraticos = sum(diferenciasPosiciones.^2, 2); % Suma de cuadrados de diferencias en x e y
ATE = sqrt(mean(erroresCuadraticos)); % Raíz cuadrada del error cuadrático medio (RMSE)

% Mostrar el ATE
disp(['ATE: ', num2str(ATE), ' metros']);

% Plot sin incertidumbre
figure;
hold on;
plot(trayectoriaReal(:,1), trayectoriaReal(:,2), 'b');
plot(trayectoriaEstimada(:,1), trayectoriaEstimada(:,2), 'r--');
hold off;
legend('Trayectoria Real', 'Trayectoria Estimada');
xlabel('Posición X (m)');
ylabel('Posición Y (m)');
title('Comparación de Trayectorias');
grid on;
saveas(gcf, 'trayectoria_sin_incertidumbre.epsc', 'epsc');

% Plot con incertidumbre
figure;
hold on;
plot(trayectoriaReal(:,1), trayectoriaReal(:,2), 'b');
plot(trayectoriaEstimada(:,1), trayectoriaEstimada(:,2), 'r--');

% Graficar las elipses de incertidumbre
for k = 1:numIteraciones
    Pk = Pk_list{k};
    [V, D] = eig(Pk(1:2, 1:2));
    t = linspace(0, 2*pi);
    chi2_critical_95 = chi2inv(0.95, 2);
    factor_escala = sqrt(chi2_critical_95);
    a = sqrt(D(1,1))*factor_escala;
    b = sqrt(D(2,2))*factor_escala;
    ellipse_x_r = a*cos(t);
    ellipse_y_r = b*sin(t);
    ellipse = [ellipse_x_r;ellipse_y_r]' * V';
    plot(ellipse(:,1)+trayectoriaEstimada(k,1), ellipse(:,2)+trayectoriaEstimada(k,2), 'r');
end

hold off;
legend('Trayectoria Real', 'Trayectoria Estimada', 'Incertidumbre Estimada');
xlabel('Posición X (m)');
ylabel('Posición Y (m)');
title('Comparación de Trayectorias con Incertidumbre');
grid on;
saveas(gcf, 'trayectoria_con_incertidumbre.epsc', 'epsc');


