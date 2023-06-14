function converge = conv(file)
shmip = load(file)
N = [shmip.md.results.TransientSolution.EffectivePressure]/(1e6);
Q = [shmip.md.results.TransientSolution.ChannelDischarge];
tt = [shmip.md.results.TransientSolution.time];

dN = N(:, 2:end) - N(:,1:end-1);
dQ = Q(:, 2:end) - Q(:, 1:end-1);
dt = tt(2:end) - tt(1:end-1);

dQdt = dQ/dt;
dNdt = dN/dt;



if abs(dQdt(end)) < 1e-4
    disp('Q is in SS')
end

