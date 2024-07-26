%Different Freqencies in ACD
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

%% Calculate Baselines (ACD w/ wf_none, dtheta_max=0.0005*2*pi
rng('default')
[d0, tvec] = test_bench(1300, 0.0005*2*pi, 130, @wf_none);

baseline_rho_bar=sum(d0.model.rho(d0.nstart:end)) / (d0.model.nsamples - d0.nstart + 1);
baseline_energy=sum(d0.trg)/d0.model.fs;

%% Run our simulation
frequency_space=linspace(70, 190, 7);
rho_bar=zeros(size(frequency_space));
energy=zeros(size(frequency_space));

for i=1:length(frequency_space)
    rng('default')

    [d1, tvec] = test_bench(1300, 0.0005*2*pi, frequency_space(i), @wf_none);
    rho_bar(i)=sum(d1.model.rho(d1.nstart:end)) / (d1.model.nsamples - d1.nstart + 1);
    energy(i)=sum(d1.trg)/d1.model.fs;

    figure
    d1.plot_osc_trg
end

% Make plots
figure
px=zeros(1,2);

px(1)=subplot(2,1,1);
plot(frequency_space, rho_bar);
yline(baseline_rho_bar, '-', 'Baseline');
xlabel('Max Frequency');
ylabel('$\rho$', 'Interpreter', 'latex');

px(2)=subplot(2,1,2 );
plot(frequency_space, energy);
yline(baseline_energy, '-', 'Baseline');
xlabel('Max Frequency');
ylabel('Energy');