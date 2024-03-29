% Useful start and end positions:
%     DESDE COCINA POR LOS PASILLOS ESTRECHOS
%         startPositionx = 8.8;
%         startPositiony = 3;
%         startOrientation = pi/2;
%         endPositionx = 20; %añadir en listaObjetivos
%         endPositiony = 13;
    
%     OBSTÁCULOS DELANTE
%         startPositionx = 5.5;
%         startPositiony = 5.5;
%         startOrientation = pi/2;
%         endPositionx = 5.5; %añadir en listaObjetivos
%         endPositiony = 16;

%     DEMOSTRADOR FINAL
%         startPositionx = 8.8;
%         startPositiony = 3;
%         startOrientation = pi/2;
%         % End positions
%         [4, 9]; % Mesa 1
%         [13.5, 17.6]; % Mesa 44
%         [13.65, 8.7]; % Mesa 69
%         [8.8, 3]; % Cocina


%%% START SETUP %%%
clear
close all
% Mapa de celdas del entorno
nombre_archivo = 'mapa_ampliado2.txt';

map = dlmread(nombre_archivo);
alturaMapa = size(map, 1);

% Costes para el A*
costs = ones(size(map));

%Punto inicial (x, y) en coordenadas del SIMULADOR
startPositionx = 8.8;
startPositiony = 3;
startOrientation = pi/2;

% Definición de la lista de puntos objetivo en coordenadas del SIMULADOR
listaObjetivos = [
    [4, 9]; % Mesa 1
    [13.5, 17.6]; % Mesa 44
    [14, 9]; % Mesa 69
    [8.8, 3]; % Cocina

    % Agregar más objetivos aquí si es necesario
];

% Convertir los objetivos a coordenadas precisas del simulador
objetivosPrecisos = round(4 * listaObjetivos);

% Índice para controlar cuál objetivo se está persiguiendo actualmente
indiceObjetivoActual = 1;

% Obtiene el primer objetivo de la lista
objetivoFinal = objetivosPrecisos(indiceObjetivoActual, :);

% Inicialización de la posición inicial y la covarianza
Xk = [startPositionx; startPositiony; startOrientation]; % Posición inicial del robot [x, y, theta]
apoloPlaceMRobot('Camarero', [startPositionx startPositiony 0], startOrientation);
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


% Parámetros para la evasión de obstáculos
distanciaEvitarObstaculo = 0.4; % Distancia para empezar a evitar obstáculos 0.3
distanciaEvitarObstaculoLateral = 0.5; % Distancia para empezar a evitar obstáculos 0.5
ultimaPosicion = [0, 0];
ultimaActualizacion = 0; % temporizador
distanciaUmbral = 0.1; % Distancia mínima esperada que el robot debe mover
tiempoUmbral = 100; % Tiempo máximo en iteraciones del bucle principal para moverse distanciaUmbral
tiempoRecuperacion = 30;
atrapado = false;
BufferFinDeEsquinaIzq = 0;
BufferFinDeEsquinaDer = 0;

%%% END SETUP %%%

%%% START BUCLE PRINCIPAL %%%
numIteraciones = 0;
k = 0;
while indiceObjetivoActual <= size(objetivosPrecisos, 1)
    k = k + 1;
    numIteraciones = numIteraciones + 1;

    %%%% START PLANIFICACIÓN %%%
    % Convertir la posición actual del robot en coordenadas del mapa
    posSegunSimulador = Xk(1:2)';
    posSegunMatriz = [4*posSegunSimulador(1), alturaMapa - 4*posSegunSimulador(2)+1]; % Invertir 'y'
    theta = Xk(3);
    
    % Nuevo objetivo final
    objetivoFinalMatriz = [objetivoFinal(1), alturaMapa - objetivoFinal(2)+1];
    %indiceTrayectoria = 1; % Indice para recorrer la trayectoria

    % Convertir coordenadas (x, y) a índices lineales para A*
    start = sub2ind(size(map), round(posSegunMatriz(2)), round(posSegunMatriz(1)));  % y from top and x from the left
    goal = sub2ind(size(map), objetivoFinalMatriz(2), objetivoFinalMatriz(1));

    % Calcular la trayectoria con A*
    if k == 1 %DEBUG
        final = a_star(logical(map), costs, start, goal);
        finalNew = final;
        a_star_plot(logical(map), costs, final); %DEBUG
        dummy=1; %DEBUG, poner breakpoint aquí para visualizar el mapa
    end %DEBUG

    final = a_star(logical(map), costs, start, goal);

    if mod(k, 500) == 0 %DEBUG, presenta mapa nuevo cada 500 iteraciones
        a_star_plot(logical(map), costs, final); %DEBUG
        dummy=1; %DEBUG
    end %DEBUG

    % Convertir los índices lineales a coordenadas (x, y)
    final_coordinates = [];
    final_flipped = flip(final);
    for i = 1:length(final_flipped)
        [y, x] = ind2sub(size(map), final_flipped(i));
        final_coordinates(i, :) = [x, alturaMapa - y+1]; % Invertir 'y' de nuevo
    end
    
    % Usar final_coordinates como la nueva trayectoriaDeseada
    % Remover la primera tupla de coordenadas de final_coordinates
    if size(final_coordinates, 1) > 1
        trayectoriaDeseada = final_coordinates(2:end, :);
    else
        % Handle the case where final_coordinates only has one element
        trayectoriaDeseada = final_coordinates;
    end
    %indiceTrayectoria = 1; % Reiniciar el índice para la nueva trayectoria

    puntoObjetivo = trayectoriaDeseada(1, :); % indiceTrayectoria = 1 siempre

    %%%% END PLANIFICACIÓN %%%

    %%%% START CONTROL REACTIVO %%%
    % Inicializar la bandera de obstáculo falso
    % Leer datos de los sensores ultrasónicos
    datosUc0 = apoloGetUltrasonicSensor('uc0'); % Sensor frontal
    datosUl1 = apoloGetUltrasonicSensor('ul1'); % Sensor izquierdo
    datosUr1 = apoloGetUltrasonicSensor('ur1'); % Sensor derecho

    % Ángulo de ajuste para los sensores laterales (en radianes)
    ajusteAngulo = 40 * (pi / 180);  % Convertir 40 grados a radianes
    
    % Calcular las nuevas orientaciones para los sensores laterales
    thetaDerecha = theta - ajusteAngulo;
    thetaIzquierda = theta + ajusteAngulo;

    % Calcular la distancia hasta la pared de la primera celda ocupada
    distanciaSensorVirtualCentro = calcularDistanciaParedCelda(map, 4*posSegunSimulador(1), alturaMapa-4*posSegunSimulador(2)+1, theta);
    distanciaSensorVirturalIzquierda = calcularDistanciaParedCelda(map, 4*posSegunSimulador(1), alturaMapa-4*posSegunSimulador(2)+1, thetaIzquierda) - 1;
    distanciaSensorVirtualDerecha = calcularDistanciaParedCelda(map, 4*posSegunSimulador(1), alturaMapa-4*posSegunSimulador(2)+1, thetaDerecha) - 1;    

    distanciaMinRealVirtualCentro = min(datosUc0, distanciaSensorVirtualCentro);
    distanciaMinRealVirtualIzquierda = min(datosUl1, distanciaSensorVirturalIzquierda);
    distanciaMinRealVirtualDerecha = min(datosUr1, distanciaSensorVirtualDerecha);

    % Verificar si hay obstáculos cercanos
    obstaculoDelante = distanciaMinRealVirtualCentro < distanciaEvitarObstaculo;
    obstaculoIzquierda = distanciaMinRealVirtualIzquierda < distanciaEvitarObstaculoLateral;
    obstaculoDerecha = distanciaMinRealVirtualDerecha < distanciaEvitarObstaculoLateral;
    
    %%% START MECANISMO DE RECUPERACIÓN DE ATRAPAMIENTO %%%

    % Calcular la distancia desde la última posición conocida
    distanciaRecorrida = norm(posSegunSimulador - ultimaPosicion);
    if distanciaRecorrida >= distanciaUmbral
        % Actualizar la última posición y reiniciar el temporizador
        ultimaPosicion = posSegunSimulador;
        ultimaActualizacion = k;
    end

    % Si el tiempo excede el umbral y la distancia no supera el umbral de movimiento,
    % entonces se considera que el robot está atrapado
    if k >= ultimaActualizacion + tiempoUmbral && distanciaRecorrida < distanciaUmbral
        % El robot no se ha movido más allá del umbral en el tiempo esperado
        atrapado = true;
        % Ejecutar procedimiento de recuperación
    end
    
    % Si atrapado, ejecutar el retroceso
    if atrapado
        disp('El robot está atrapado y se está ejecutando el procedimiento de recuperación.');
        BufferFinDeEsquinaIzq = 0;
        BufferFinDeEsquinaDer = 0;
        comandoVelocidadLinear = -0.1;
        comandoVelocidadAngular = 0;
        tiempoRecuperacion = tiempoRecuperacion - 1;
        if tiempoRecuperacion == 0
            tiempoRecuperacion = 30;
            atrapado = 0;
        end
    end
    
    %%% END MECANISMO DE RECUPERACIÓN DE ATRAPAMIENTO %%%
  
    % Si hay un obstáculo, ajustar la velocidad y la dirección
    if (obstaculoDelante || obstaculoIzquierda || obstaculoDerecha) && ~atrapado
        BufferFinDeEsquinaIzq = 0;
        BufferFinDeEsquinaDer = 0;
        
        if obstaculoDelante

            % Obstáculo delante e izquierda
            if obstaculoIzquierda && ~obstaculoDerecha
                % Ajustar la dirección a la derecha
                comandoVelocidadLinear = 0.02;
                comandoVelocidadAngular = -Kp_angular * 2;
            end
            
            % Obstáculo delante y derecha
            if obstaculoDerecha && ~obstaculoIzquierda
                % Ajustar la dirección a la izquierda
                comandoVelocidadLinear = 0.02;
                comandoVelocidadAngular = Kp_angular * 2;                
            end
            
            % Sólo obstaculo delante. Es poco ancho.
            if ~obstaculoIzquierda && ~obstaculoDerecha
                % Parar y girar
                comandoVelocidadLinear = 0; % Detenerse o reducir velocidad
                if distanciaMinRealVirtualIzquierda < distanciaMinRealVirtualDerecha
                    comandoVelocidadAngular = -Kp_angular; % Girar a la derecha
                elseif distanciaMinRealVirtualIzquierda > distanciaMinRealVirtualDerecha
                    comandoVelocidadAngular = Kp_angular; % Girar a la izquierda
                else
                    comandoVelocidadAngular = -Kp_angular; % Girar a la derecha or defecto
                end
            end

            % Obstáculo delante, derecha e izquierda
            if obstaculoIzquierda && obstaculoDerecha
                % Reducir la velocidad lineal y girar
                comandoVelocidadLinear = 0.01; % Detenerse o reducir velocidad
                if distanciaMinRealVirtualIzquierda < distanciaMinRealVirtualDerecha
                    comandoVelocidadAngular = -Kp_angular; % Girar a la derecha
                elseif distanciaMinRealVirtualIzquierda > distanciaMinRealVirtualDerecha
                    comandoVelocidadAngular = Kp_angular; % Girar a la izquierda
                else
                    comandoVelocidadAngular = -Kp_angular; % Girar a la derecha or defecto
                end
            end

        % Obstáculo izquierda y derecha, hay un hueco
        elseif obstaculoIzquierda && obstaculoDerecha
            % Error entre la posición actual estimada y el punto objetivo
            error = puntoObjetivo - 4*(Xk(1:2)');
            
            % Controlador proporcional
            Kp = 1; % Ganancia proporcional
            comandoVelocidadLinear = Kp * norm(error);
            
            % Calcular la velocidad angular necesaria para girar hacia el punto objetivo
            anguloObjetivo = atan2(error(2), error(1));
            errorAngular = anguloObjetivo - Xk(3);
            % Asegurarse de que el error angular esté en el rango [-pi, pi]
            errorAngular = mod(errorAngular + pi, 2 * pi) - pi;
            
            Kp_angular = 0.5; % Ganancia proporcional para la velocidad angular
            comandoVelocidadAngular = Kp_angular * errorAngular;

        % Sólo obstaculo izquierda
        elseif obstaculoIzquierda 
            % Ajustar la dirección ligeramente a la derecha
            comandoVelocidadLinear = 0.05;
            comandoVelocidadAngular = -Kp_angular * 1;

            % Buffer para evitar el problema de fin de esquina
            BufferFinDeEsquinaIzq = 40;

        % Sólo obstaculo derecha    
        elseif obstaculoDerecha 
            % Ajustar la dirección ligeramente a la izquierda
            comandoVelocidadLinear = 0.05;
            comandoVelocidadAngular = Kp_angular * 1;

            % Buffer para evitar el problema de fin de esquina
            BufferFinDeEsquinaDer = 40;

        end


    % No hay obstaculos detectados 
    elseif ~atrapado

        % Error entre la posición actual estimada y el punto objetivo
        error = puntoObjetivo - 4*(Xk(1:2)');
        
        % Controlador proporcional
        Kp = 1; % Ganancia proporcional
        comandoVelocidadLinear = Kp * norm(error);
        
        % Calcular la velocidad angular necesaria para girar hacia el punto objetivo
        anguloObjetivo = atan2(error(2), error(1));
        errorAngular = anguloObjetivo - Xk(3);
        % Asegurarse de que el error angular esté en el rango [-pi, pi]
        errorAngular = mod(errorAngular + pi, 2 * pi) - pi;
        
        Kp_angular = 0.5; % Ganancia proporcional para la velocidad angular
        comandoVelocidadAngular = Kp_angular * errorAngular;
        
        %%% START BUFFER FIN DE ESQUINA %%%
        if BufferFinDeEsquinaIzq
            BufferFinDeEsquinaIzq = BufferFinDeEsquinaIzq - 1;
            if BufferFinDeEsquinaIzq == 0 %DEBUG
                %a_star_plot(logical(map), costs, final); %DEBUG
            end %DEBUG
            if comandoVelocidadAngular > 0
                comandoVelocidadAngular = 0;
            end
        end
    
        if BufferFinDeEsquinaDer
            BufferFinDeEsquinaDer = BufferFinDeEsquinaDer - 1;
            if BufferFinDeEsquinaDer == 0 %DEBUG
                %a_star_plot(logical(map), costs, final); %DEBUG
            end %DEBUG
            if comandoVelocidadAngular < 0
                comandoVelocidadAngular = 0;
            end
        end
        %%% END BUFFER FIN DE ESQUINA %%%
    else
        dummy = 1;
    end

    

    % Limitar la velocidad máxima lineal
    velocidadMaxLineal = 0.1; % Velocidad máxima lineal
    comandoVelocidadLinear = min(comandoVelocidadLinear, velocidadMaxLineal);

    % Limitar la velocidad máxima angular
    velocidadMaxAngular = pi/4; % Velocidad máxima angular, por ejemplo, pi/4 rad/s


    comandoVelocidadAngular = max(min(comandoVelocidadAngular, velocidadMaxAngular), -velocidadMaxAngular);
    timestep = 0.1;

    % Ejecutar el comando de movimiento
    apoloMoveMRobot('Camarero', [comandoVelocidadLinear, comandoVelocidadAngular], timestep);
    apoloUpdate();

    %%%% END CONTROL REACTIVO %%%
    

    %%%% START LOCALIZACIÓN %%%
    % Obtención de la posición real del robot con apoloGetLocation
    posSegunSimulador = apoloGetLocation('Camarero'); % Esto devuelve un vector de 6 elementos
    
    % Seleccionamos solo los valores de interés: x, y y yaw
    posReal = [posSegunSimulador(1); posSegunSimulador(2); posSegunSimulador(6)]; % Ahora posReal tendrá solo x, y, yaw

    % Obtención de las observaciones de las balizas visibles
    laserLandmarks = apoloGetLaserLandMarks('LMS100');

    % Obtención de Uk con apoloGetOdometry
    odometry = apoloGetOdometry('Camarero');

    % Verificar si odometry(1) es negativo
    if odometry(1) < 0
        % Si es negativo, calcular Uk con signo negativo en el desplazamiento
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
        %disp(mahalanobisDistances)
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
        %%%% END LOCALIZACIÓN %%%
    end

    % Verificar si el robot ha alcanzado el objetivo actual
    if isequal(round(4 * (Xk(1:2)')), objetivoFinal)
        % Avanzar al siguiente objetivo
        indiceObjetivoActual = indiceObjetivoActual + 1;
        if indiceObjetivoActual <= size(objetivosPrecisos, 1)
            objetivoFinal = objetivosPrecisos(indiceObjetivoActual, :);
            % ... [Código para recalcular la trayectoria al nuevo objetivo] ...
        else
            % Todos los objetivos alcanzados
            disp("Todos los objetivos alcanzados");
            break;
        end
    end
end
%%% END BUCLE PRINCIPAL %%%

%%% START GRAFICAS %%%
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

%%% END GRAFICAS %%%
