function qgrn_draw_network(genes, theta, threshold, filename)
%QGRN_DRAW_NETWORK  Plot the fitted gene regulatory network.
%
%   QGRN_DRAW_NETWORK(GENES, THETA, THRESHOLD, FILENAME) draws one node per
%   gene on a circle and one edge per surviving off-diagonal angle.
%
%   Pruning rule from the paper: an edge is kept only if
%       |theta| > THRESHOLD * pi/180        (THRESHOLD = 1 degree by default)
%   Edge width is 30|theta|/pi, green for theta >= 0 (activation) and red for
%   theta < 0 (repression) -- the same encoding as the Python draw_network.

    if nargin < 3 || isempty(threshold), threshold = 1; end
    if nargin < 4, filename = ''; end

    n = numel(genes);
    cut = threshold * pi / 180;

    s = []; t = []; w = []; sgn = [];
    for a = 1:n
        for b = a+1:n                                  % theta is symmetric
            if abs(theta(a,b)) > cut
                s(end+1) = a;  t(end+1) = b;                      %#ok<AGROW>
                w(end+1) = 30 * abs(theta(a,b)) / pi;             %#ok<AGROW>
                sgn(end+1) = sign(theta(a,b));                    %#ok<AGROW>
            end
        end
    end

    ang = (0:n-1) * 2*pi/n;
    fig = figure('Color', 'w', 'Position', [100 100 560 560]);
    if isempty(s)
        text(0.5, 0.5, 'no edges above threshold', 'HorizontalAlignment', 'center');
        axis off;
    else
        G = graph(s, t, w, cellstr(genes));
        h = plot(G, 'XData', cos(ang), 'YData', sin(ang), ...
                 'LineWidth', max(w/max(w)*6, 0.5), ...
                 'NodeColor', [0.855 0.890 0.949], ...
                 'MarkerSize', 26, 'NodeFontSize', 11);
        ec = zeros(numel(sgn), 3);
        ec(sgn >= 0, :) = repmat([0.439 0.678 0.278], sum(sgn >= 0), 1);  % green
        ec(sgn <  0, :) = repmat([1.000 0.000 0.000], sum(sgn <  0), 1);  % red
        h.EdgeColor = ec;
        axis equal off;
        title(sprintf('QuantumGRN  (%d edges, |\\theta| > %.1f\\circ)', ...
                      numel(s), threshold));
    end

    if ~isempty(filename)
        exportgraphics(fig, filename, 'ContentType', 'vector');
        fprintf('  wrote %s\n', filename);
    end
end
