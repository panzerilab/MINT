function [B_tem, H, B_spa, err] = sNM3F(R, P, L, opts)
% Space by time

if nargin < 4
    opts.max_iter = 5000;
    opts.tolerance = 0.01;
end

[N, T, S] = size(R);

% Step 1: Initialize B_tem, H and B_spa with random entries
B_tem = rand(T, P);
H = rand(P, L, S);
B_spa = rand(L, N);

err_all=NaN(1,opts.max_iter);

for i = 1:(opts.max_iter)
   
    [B_tem, B_spa, H] = normalize(B_tem, B_spa, H, S, P, L);
    % Step 2: update B_spa
    H_reshaped = reshape(H, [P, L*S]);
    G = B_tem*H_reshaped;
    G_reshaped = reshape(G, [T*S, L]);
    R_reshaped = reshape(R, [T*S, N]);

    numerator = G_reshaped' * R_reshaped;
    denominator = (G_reshaped' * G_reshaped) * B_spa;
    B_spa = B_spa .* (numerator ./ (denominator + eps(numerator)));

    % Step 3: update B_tem
    H_reshaped = reshape(H, [P*S,L]);
    V = H_reshaped*B_spa;
    V_reshaped = reshape(V, [P, N*S]);
    R_reshaped = reshape(R, [T, N*S]);

    numerator = R_reshaped * V_reshaped';
    denominator = B_tem * (V_reshaped * V_reshaped');
    B_tem = B_tem .* (numerator ./ (denominator + eps(numerator)));

    % Step 4: update H
    for s = 1:S
        R_s = squeeze(R(:, :, s))';
        H_s = squeeze(H(:, :, s));

        numerator = B_tem' * R_s * B_spa';
        denominator = (B_tem' * B_tem) * H_s * (B_spa * B_spa');
        H(:, :, s) = H_s .* (numerator ./ (denominator + eps(numerator)));
    end

    % Step 5: calculate error
    err = calculateErrorSum(R, B_tem, H, B_spa, S);
    err_all(i) = err;

    if i>2           
        if err_all(i)>err_all(i-1)
            break;
        end
    end

    if err < opts.tolerance
        [B_tem, B_spa, H] = normalize(B_tem, B_spa, H, S, P, L);
        return
    end
end
end

function err = calculateErrorSum(R, B_tem, H, B_spa, S)
err = 0;
for s = 1:S
    H_s = squeeze(H(:, :, s));
    reconstruction_s = (B_tem * H_s * B_spa)';
    err = err + norm(R(:, :, s) - reconstruction_s, 'fro')^2;
end
end

function [B_tem, B_spa, H] = normalize(B_tem, B_spa, H, S, P, L)
alpha = sqrt(sum(B_tem.^2, 1));
beta = sqrt(sum(B_spa.^2, 2));

B_tem = B_tem ./ alpha;
B_spa = B_spa ./ beta;

alpha_inv = 1 ./ alpha;
beta_inv = 1 ./ beta;
for s = 1:S
    for p = 1:P
        for l = 1:L
            H(p, l, s) = H(p, l, s) .* alpha(p)*beta(l);
        end
    end
end
end
