%Simple shape testing setup
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

%% Simulate

rng(45)
frequency=5;
[d1, tvec] = test_bench(400, 0.0005*2*pi*130/frequency, frequency, @wf_none);



%% Plot data

figure
%Plot symptom signal together with average stimulation
%triggers across stimulating electrode.
d1.plot_osc_energy;

%figure
%Plot triggers across stimulating electrodes.
%d1.plot_trg_l;

%figure
%Plot electrode channel activity.
% d1.plot_channels;

%% Show system

%figure
% plot_system_spherical;

