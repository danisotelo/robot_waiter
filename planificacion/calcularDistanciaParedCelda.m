% Función para calcular la distancia hasta el borde de la primera celda ocupada
function distancia = calcularDistanciaParedCelda(map, posX, posY, theta)
    dx = cos(theta);
    dy = -sin(theta);  % '-sin' debido a la inversión de 'y'
    distancia = inf;  % Inicializar con un valor alto

    % Iterar hasta encontrar una celda ocupada
    for d = 0:0.01:inf  % Asumiendo que cada celda tiene una unidad de longitud
        x = posX + d * dx;
        y = posY + d * dy;

        % Encontrar la celda actual
        celdaX = floor(x);
        celdaY = floor(y);

        % Comprobar si estamos fuera de los límites del mapa
        if celdaX <= 0 || celdaX > size(map, 2) || celdaY <= 0 || celdaY > size(map, 1)
            break;
        end

        % Comprobar si la celda está ocupada
        if map(celdaY, celdaX) == 0  % COLUMN; ROW OF MATRIX. NO PUEDEN SER 0, a partir del 1.
            distancia = sqrt((x - posX)^2 + (y - posY)^2);  % Calcular la distancia euclidiana
            break;
        end
    end
end
