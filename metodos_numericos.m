%{
Proyecto final de Métodos Numericos, semestre Febrero - Junio 2020 

Profesor: Jorge Islas
Integrantes:
    * Diego de Jesús Ramírez Rodríguez 
    * Alejandra Velázquez Mendoza

repo: https://github.com/DiegoJRR/proyecto-MN

Escrito durante la 'Plaga de los 20'

Última revisión: 08/22/20
%}


% Introducción
disp("Este es el proyecto de métodos numericos de Diego Ramírez y Alejandra Velázquez")
disp_with_space("El proyecto contiene métodos para reslvoer sistemas lineales de ecuaciones y funciones de una variable")

disp("Para comenzar, ingrese una matriz válida para los métodos de solución de sistemas lineaes")
disp("Para el mejor resultado posible, ingrese una matriz diagonalmente dominante")
response = input("Si quiere utilizar la matriz default (y/n)", 's');

if response == "y"
    A = [3 2 1 -2 -6; 1 5 2 1 -1; 2 1 6 1 -5 ; 3 1 -2 9 -1 ; 2 1 3 4 19]
    [n, m] = size(A);
else
    disp("Debe estar en notación de MATLAB: ´[fila_1 ; fila_2 ; ... ; fila_n]´ donde las filas se escriben de forma: 'a_1 a_2 a_3' con los valores separados por espacios")

    A = input("Ingrese la matriz cuadrada con el formato indicado: ")
    [n, m] = size(A);

    while isempty(A) || n ~= m
        disp("Para comenzar, ingrese una matriz válida para los métodos de solución de sistemas lineaes")
        disp("Debe estar en notación de MATLAB: ´[fila_1 ; fila_2 ; ... ; fila_n]´ donde las filas se escriben de forma: 'a_1 a_2 a_3' con los valores separados por espacios")

        A = input("Ingrese la matriz con el formato indicado: ")
        [n, m] = size(A);
    end

    disp(" ")
end

disp("Ingrese ahora un vector que quiere utilizar como la solución del sistema A*x")
disp("El formato para los vectores es: '[a_1 a_2 ... a_n]', asegurese que tiene el mismo número de filas que la matriz cuadrada")

b = input("Ingrese el vector con el formato indicado: ")
[n_b, m_b] = size(b);

while isempty(b) || m_b ~= n || n_b ~= 1
    disp("Ingrese ahora un vector que quiere utilizar como la solución del sistema A*x")
    disp("El formato para los vectores es: '[a_1 a_2 ... a_n]', asegurese que tiene el mismo número de filas que la matriz cuadrada")

    b = input("Ingrese el vector con el formato indicado: ")
    [n_b, m_b] = size(b);
end



disp(" ")

% Sección con métodos matriciales
disp_with_space("Los siguientes son los métodos matriciales para resolver sistemas de ecuaciones, utilizando la matriz dada arriba")

disp("El primero es el método de Gauss, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

input("ENTER para continuar")

x = solve_with_gauss(A, b);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)


%  Solución con método de Cramer
disp("El siguiente es el método de Cramer, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

input("ENTER para continuar")

x = solve_with_cramer(A, b);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)


%  Solución con método de Gauss Jordan
disp("El siguiente es el método de Gauss Jordan, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

input("ENTER para continuar")

x = solve_with_gauss_jordan(A, b);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)



%  Solución con método de Matriz Inversa
disp("El siguiente es el método de Matriz Inversa, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

input("ENTER para continuar")

x = solve_with_inverse(A, b);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)


%  Solución con método de LU
disp("El siguiente es el método de LU, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

input("ENTER para continuar")

x = solve_with_LU(A, b);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)


%  Solución con método de Guass Seidel
disp("El siguiente es el método de Gauss Seidel, utilizando la matriz")
disp(A)
disp("Y el vector b")
disp(b)

threshold = input("Ingrese el porcentaje máximo de precisión que requiere: ")/100;

input("ENTER para continuar")

x = solve_with_gauss_seidel(A, b, threshold);

disp("Obtenemos el resultado: ")
disp_with_space(x)

disp("Y multiplicando A*x vemos que obtenemos b de nuevo")
disp_with_space(A*x)

% Sección de métodos numericos para funciones univariables

disp("Los siguientes métodos numericos encuentran las raices de funciones univariables constantes")
disp("Ingrese una función matemática con la notación de MATLAB")

syms x
f(x) = str2sym(input("Ingrese una función univarible de x: ", 's'))
threshold = input("Ingrese un valor máximo del error deseado (porcentual): ")/100;

disp(" ")
disp(" ")

% Solución con método de Bisección
disp_with_space("Solución con el método de Bisección")

x_l = input("Ingrese el valor de x_l para el método: ")
x_u = input("Ingrese el valor de x_u para el método: ");

x = solve_with_bisection(f, x_l, x_u, threshold);

disp("Obtenemos entonces la solución")
disp(x)
disp("Que corresponde a un valor de f")
disp_with_space(f(x))


% Solución con método de Newton - Raphson
disp_with_space("Solución con el método de Newton - Raphson")

x_0 = input("Ingrese el valor de x_0 para el método: ");

x = solve_with_newton_rapshon(f, x_0, threshold);

disp("Obtenemos entonces la solución")
disp(x)
disp("Que corresponde a un valor de f")
disp_with_space(f(x))


disp("Gracias por utilizar este programa. Este archivo incluye todas las funciones utilizadas de los métodos, y algunas otras funciones de apoyo, y son compatibles independientemente con objetos de MATLAB.")


function disp_with_space(text)
    disp(text)
    input("ENTER para continuar")
    disp(" ")
    disp(" ")
end

function x = solve_with_newton_rapshon(f, x_0, threshold)
    x_approx = zeros(2, 1);
    x_approx(1) = x_0;
    df = diff(f);
    
    error = 1;
    
    while error > threshold
        temp = x_approx(2);
        x_approx(2) = x_approx(1) - f(x_approx(1))/df(x_approx(1));
        
        x_approx(1) = temp;
        
        error = abs((x_approx(2) - x_approx(1))/x_approx(2));
    end
    
    x = x_approx(2);
end

function x = solve_with_bisection(f, x_l, x_u, threshold)  
    if f(x_l)*f(x_u) > 0
        disp("Los valores iniciales deben tener signos distintos")
        x = false
        return
    end
    
    % Define a vector with x_r(i) being the old approximation
    x_r = zeros(2, 1);
    
    error = 1;
    
    while error > threshold
        x_r(1) = x_r(2);
        x_r(2) = (x_l + x_u)/2;
        
        error = abs((x_r(2) - x_r(1))/x_r(2));
        product = f(x_l)*f(x_r(2));
        
        if product < 0
            x_u = x_r(2);
            continue
        elseif product > 0
            x_l = x_r(2);
            continue
        else
            x = x_r(2);
            return
        end
        
        
    end
    
    x = x_r(2);
end

function sol = solve_with_gauss_seidel(A, b, threshold)
    [n, ~] = size(A);
    
    if ~check_diagonally_dominant(A)
        disp("No se puede resolver con este metodo porque la matriz no es diagonalmente dominante")
        sol = false
        return
    end
    
    x = zeros(n, 1);
    errors = zeros(n, 1);
    
    for i = 1:n
        errors(i) = 100;
    end
    
    while sum(errors < threshold) ~= n
        for i = 1:n
            aux = x;
            x_last = aux(i);
            aux(i) = 0;
            
            % Update x_i value 
            x(i) = (b(i) - dot(aux, A(i, 1:n)))/A(i, i);
           
            % Update errors_i value
            errors(i) = abs((x(i) - x_last)/x(i))*100;
        end
        
    end
    
    sol = x;
end

function x = check_diagonally_dominant(A)
    [n, ~] = size(A);
    
    for row = 1:n
       if abs(A(row, row)) >= sum_abs(A(row, :)) - abs(A(row, row))
           continue
       else
           x = false;
           return
       end
    end
    
    x = true;
end

function x = sum_abs(v)
    s = 0
    for i = 1:size(v)
        s = s + abs(v(i));
    end
    
    x = s;
end

function x = solve_with_cramer(A, b)
    [n, ~] = size(A);
    det_A = determinant(A);
    
    if det_A == false || det_A == 0
        x = false;
        return
    end
    
    x = zeros(n, 1);
    
    for i = 1:n
       aux = A;
       aux(:, i) = b;
       
       x(i) = determinant(aux)/det_A;
    end
end

function x = solve_with_LU(A, b)
    n = length(b);
    x = zeros(n, 1);
    y = zeros(n, 1);

    for i = 1:1:n
        for j = 1:1:(i - 1)
            a = A(i,j);
            for k = 1:1:(j - 1)
                a = a - A(i,k)*A(k,j);
            end
            A(i,j) = a/A(j,j);
        end
        for j = i:1:n
            a = A(i,j);
            for k = 1:1:(i - 1)
                a = a - A(i,k)*A(k,j);
            end
            A(i,j) = a;
        end
    end

    for i = 1:1:n
        a = 0;
        for k = 1:1:i
            a = a + A(i,k)*y(k);
        end
        y(i) = b(i) - a;
    end

    for i = n:(-1):1
        a = 0;
        for k = (i + 1):1:n
            a = a + A(i,k)*x(k);
        end
        x(i) = (y(i) - a)/A(i, i);
    end    
end

function x = solve_with_gauss_jordan(A, b)
    b = b';
    
    if determinant(A) == false
        x = false;
        return
    end
    
    aug = augment(A, b);
    U = make_upper(aug);
    
    ident_aug = make_ident(U);
    ident_aug = reduce_ref(ident_aug);
    
    x = ident_aug(:, end);
end

function x = solve_with_gauss(A, b)
    b = b';
    
    if determinant(A) == false
        x = false;
        return
    end
    
    aug = augment(A, b);
    A = make_upper(aug);
    
    x = solve_ref(A);
end

function x = solve_with_inverse(A, b)
    if determinant(A) == false
        x = false;
        return
    end
    
    [n, ~] = size(A);
    
    aug = augment(A, eye(n));
    U = make_upper(aug);
    ident_aug = make_ident(U);
    
    ident_aug = reduce_ref(ident_aug);
    inverse = ident_aug(:, n+1:n+n);
    
    x = inverse*b';
end

function aug = augment(A, b)
   A = [A, b];
 
   aug = A;
end

function sol = make_ident(A)
    % Esta función toma como input una matriz aumentada con el lado
    % izquierdo en forma superior, y lo convierte a la identidad
    [n, ~] = size(A);
    
    for row_reverse = n:-1:2
       for row = 1:row_reverse - 1
           A(row, :) = A(row, :) - (A(row, row_reverse)/A(row_reverse, row_reverse))*A(row_reverse, :);
       end
    end
    
    sol = A;
end

function sol = reduce_ref(A)
    [n, ~] = size(A);
    
   for row = 1:n       
       A(row, :) = A(row, :)/A(row, row);
   end
    
    sol = A;
end

function det = determinant(A)
    % Calcular la determinante de una matriz A
    [n, m] = size(A);
    
    if n ~= m
        det = false;
        return
    end
    
    A = make_upper(A);

    
    if A == false
        disp('La matriz no tiene solución única porque su determinante es 0');
        det = false;
        return
    end
    
    determinant_value = 1;
    
    for row = 1:n
        for col = 1:m
            if row == col
                determinant_value = determinant_value*A(row, col);
            end
        end
    end
    
    det = determinant_value;
    return
end

function sol = make_upper(A)
    % Convierte a forma fila-escalon una matrix
    % m_limit es la ultima columna antes de que comience la parte aumentada
   [n, ~] = size(A);
   
   % Asegurar que la fila 1 tiene valor diferente de 0
   for row = 1:n
       for col = 1:n
           if row == col
               
              if A == 0
                 sol = false;
                 return
              end
              
              A = make_pivot_column(A, row, col);
           end
       end
   end
   
   sol = A;
end

function sol = solve_ref(A)
    [n, m] = size(A);
    
    if n ~= m-1
        sol = false;
        return
    end
    
    b = A(:, m);
    
    x = zeros(n, 1);
    x(n) = b(n)/(A(n, n));
    
    for i = n-1:-1:1
        x(i) = (b(i) - dot(A(i, i:n), x(i:n)))/(A(i, i));
    end
    
    sol = x;
end

function x = make_pivot_column(A, i, j)
    A = ensure_pivot(A, i, j);
    
    if A == false
        x = false;
        return
    end
    
    [n, ~] = size(A);
        
    for row = i+1:n
       A(row, :) = A(row, :) - (A(row, j)/A(i, j))*A(i, :);
    end
    
    x = A;
    
end

function x = ensure_pivot(A, i, j)    
    if A(i, j) ~= 0
        x = A;
    else
        [n, ~] = size(A);
        
        for row = i+1:n
            if A(row, j) ~= 0
                x = exchange_rows(A, j, row);
                return
            end
        end
        
        x = false;
    end
end

function x = exchange_rows(A, row1, row2)
    temp = A(row2, :);
    A(row2, :) = A(row1, :);
    A(row1, :) = temp;
    
    x = A;
end

