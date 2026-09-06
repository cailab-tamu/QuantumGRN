function G = gate_defs()
%GATE_DEFS  Handles for single-qubit rotation gates.
    G.Ry = @(t) [cos(t/2) -sin(t/2); sin(t/2) cos(t/2)];
    G.Rx = @(t) [cos(t/2) -1i*sin(t/2); -1i*sin(t/2) cos(t/2)];
end
