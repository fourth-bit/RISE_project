%Oscilation Shape in ACD
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

d1, tvec = test_bench(1300, 0.0005*2*pi, 130, @wf_sinusoid);

% Make plots
figure
d1.plot_osc_rho