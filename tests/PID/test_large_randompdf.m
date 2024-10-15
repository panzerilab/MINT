% test_large_randompdf.m -- part of BROJA_2PID (https://github.com/dot-at/BROJA_2PID/)
% Usage: test_large_randompdf(x, y, z, iter)
% Where:   x    is the size of the range of X;
%          y    is the size of the range of Y;
%          z    is the size of the range of Z;
%          iter is the number of iterations (defaults to 250).

function test_large_randompdf(nX, nY, nZ, maxiter)
    if nargin < 3 || nargin > 4
        fprintf('Usage: test_large_randompdf(x, y, z, iter)\n');
        fprintf('Where:   x    is the size of the range of X;\n');
        fprintf('         y    is the size of the range of Y;\n');
        fprintf('         z    is the size of the range of Z;\n');
        fprintf('         iter is the number of iterations (defaults to 250).\n');
        return;
    end

    if nargin < 4
        maxiter = 250;  % Default number of iterations
    end

    if min([nX, nY, nZ]) < 2
        error('All sizes of ranges must be at least 2.');
    end

    if maxiter < 1
        error('# iterations must be >= 1.');
    end

    tStart = tic;  % Start overall timer
    Ti = zeros(1,maxiter)
    for iter = 1:maxiter
        fprintf('Random PDFs   with |X| = %d |Y| = %d |Z| = %d\n', nX, nY, nZ);
        fprintf('______________________________________________________________________\n');
        fprintf('Create pdf #%d\n', iter);

        % Create a random pdf
        pdf = zeros([nX, nY, nZ]);
        pts = sort(rand(1, nX * nY * nZ - 1),"descend"); % Generate random points
        pts = [pts, 0];  % Add a zero to the end of sorted list
        val = 1.0;

        % Fill the pdf map
        idx = 1;
        for x = 1:nX
            for y = 1:nY
                for z = 1:nZ
                    newval = pts(idx);
                    pdf(x, y, z) = val - newval;
                    val = newval;
                    idx = idx + 1;
                end
            end
        end

        fprintf('Run BROJA_2PID.pid().\n');
        % tic  % Start iteration timer
        pid_result = pidBROJA(pdf);  % Call the PID function (output = 0)
        % Ti(iter) = toc;  % Stop iteration timer

        fprintf('Partial information decomposition: \n');
        disp(pid_result);
        % fprintf('Time: %.6f secs\n', Ti(iter));
        fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
    end
    % toc_net = sum(Ti);
    toc_total = toc(tStart);  % Stop overall timer
    fprintf('**********************************************************************\n');
    fprintf('Average time: %.6f secs\n', (toc_total) / maxiter);
    fprintf('**********************************************************************\n');
    % fprintf('Average time net: %.6f secs\n', (toc_net) / maxiter);
end
