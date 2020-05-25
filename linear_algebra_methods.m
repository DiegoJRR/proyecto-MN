A = rand(3);
b = [4, 8, 12];

solution = solve_with_gauss(A, b)

A*solution

function sol = solve_with_gauss(A, b)
    b = b';
    
    aug = augment(A, b);
    [n_b, m_b] = size(b);
    A = ref(aug, m_b);
    
    sol = solve_ref(A);
end

function aug = augment(A, b)
   A = [A, b];
 
   aug = A;
end


function sol = reduce_ref(A)
    [n, m] = size(A);
    
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
        det = false
    end
    
    A = ref(A, 0)
    
    determinant_value = 1
    
    for row = 1:n
        for col = 1:m
            if n == m
                determinant_value = determinant_value*A(row, col);
            end
        end
    end
    det = A
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
              A = make_pivot_column(A, row, col);
              
              if A == false
                 sol = false;
              end
           end
       end
   end
   
   sol = A;
end

function sol = solve_ref(A)
    [n, m] = size(A);
    
    if n ~= m-1
        sol = false;
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
    end
    
    [n, m] = size(A);
        
    for row = i+1:n
       A(row, :) = A(row, :) - (A(row, j)/A(i, j))*A(i, :);
    end
    
    x = A;
    
end

function x = ensure_pivot(A, i, j)
    if A(i, j) ~= 0
        x = A;
    else
        [n, m] = size(A);
        
        for row = i+1:n
            if A(row, j) ~= 0
                x = exchange_rows(A, j, row);
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

