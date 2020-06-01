A = rand(3);
b = [4, 500, 8];
solution = solve_with_gauss(A, b);
A*solution

sol = solve_with_inverse(A, b);
A*sol


function sol = solve_with_gauss(A, b)
    b = b';
    
    if determinant(A) == false
        sol = false;
        return
    end
    
    aug = augment(A, b);
    A = make_upper(aug);
    
    sol = solve_ref(A);
end

function sol = solve_with_inverse(A, b)
    if determinant(A) == false
        sol = false;
        return
    end
    
    [n, ~] = size(A);
    
    aug = augment(A, eye(n));
    U = make_upper(aug);
    ident_aug = make_ident(U);
    
    ident_aug = reduce_ref(ident_aug);
    inverse = ident_aug(:, n+1:n+n);
    
    sol = inverse*b';
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
       for row = 1:row_reverse-1
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
   [n, m] = size(A);
   
   m = n;
   % Asegurar que la fila 1 tiene valor diferente de 0
   for row = 1:n
       for col = 1:m
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

