A = rand(3);
b = [4, 8, 12];
solution = solve_with_gauss(A, b)
A*solution

function sol = solve_with_gauss(A, b)
    b = b';
    
    if determinant(A) == false
        sol = false;
        return
    end
    
    aug = augment(A, b);
    [~, m_b] = size(b);
    A = ref(aug, m_b);
    
    sol = solve_ref(A);
end

function aug = augment(A, b)
   A = [A, b];
 
   aug = A;
end


function sol = reduce_ref(A)
    [n, ~] = size(A);
    
   for row = 1:n
       row_n = A(row, :);
       m_pivot = find(row_n, 1, 'first');
       
       A(row, :) = A(row, :)/A(row, m_pivot);
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
    
    A = ref(A, 0);

    
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

function sol = ref(A, b_cols)
    % Convierte a forma fila-escalon una matrix
    % m_limit es la ultima columna antes de que comience la parte aumentada
   [n, m] = size(A);
   
   m = m - b_cols;
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

