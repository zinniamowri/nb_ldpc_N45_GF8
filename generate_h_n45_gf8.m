clc; clear;

% Parameters
N  = 45;   % number of columns
dv = 2;    % column weight
dc = 3;    % row weight
q  = 8;    % GF(8)

M = N*dv/dc;   % number of rows

if mod(M,1) ~= 0
    error('Invalid parameters: N*dv must be divisible by dc');
end

M = round(M);

% -------------------------------------------------
% Step 1: Generate binary support matrix Hb
% -------------------------------------------------
Hb = zeros(M,N);

success = false;
max_trials = 1000;

for trial = 1:max_trials
    
    Hb = zeros(M,N);
    
    % Create sockets
    vn_sockets = repelem(1:N, dv);   % each variable repeated dv times
    cn_sockets = repelem(1:M, dc);   % each check repeated dc times
    
    % Random permutation of check sockets
    cn_sockets = cn_sockets(randperm(length(cn_sockets)));
    
    valid = true;
    
    for k = 1:length(vn_sockets)
        col = vn_sockets(k);
        row = cn_sockets(k);
        
        % Avoid parallel edge in same row/column position
        if Hb(row,col) == 1
            valid = false;
            break;
        end
        
        Hb(row,col) = 1;
    end
    
    % Check row/column weights
    if valid && all(sum(Hb,1)==dv) && all(sum(Hb,2)==dc)
        success = true;
        break;
    end
end

if ~success
    error('Could not generate a valid binary support matrix.');
end

disp('Binary support matrix Hb generated successfully.');
disp(['Size of Hb = ', num2str(M), ' x ', num2str(N)]);

% -------------------------------------------------
% Step 2: Convert to non-binary GF(8) matrix
% -------------------------------------------------
% Primitive polynomial example for GF(8): x^3 + x + 1


Hgf = gf(zeros(M,N), 3);

for i = 1:M
    for j = 1:N
        if Hb(i,j) == 1
            Hgf(i,j) = gf(randi([1,q-1]), 3);
        end
    end
end

disp('Non-binary parity-check matrix H over GF(8) generated.');

h = Hgf.x;   % extract integer representation (0 to 7)
save('nbldpc_45_2_3_gf8.mat', 'h');