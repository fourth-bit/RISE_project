%Ensembling Everything into One Script
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

%% Setup
%waveform_names={'Square', 'Expontial Rise', 'Exponential Decay', 'Sinusoidal', 'Triangular'};
%waveforms={@wf_none, @wf_exp, @wf_inv_exp, @wf_sinusoid, @wf_triangle};
waveform_names={'Square' 'Sinusodal'};
waveforms={@wf_none, @wf_sinusoid};

dtheta_maxes=linspace(0, 0.0005, 5);
freqs=linspace(90, 150, 4);

runs=8;
rhos=zeros(length(waveforms), length(dtheta_maxes), length(freqs), runs);
energies=zeros(length(waveforms), length(dtheta_maxes), length(freqs), runs);


%% Run Ensemble
tic

% For reproduction
rng('default')

i=0;
max_i=length(waveforms) * length(dtheta_maxes) * length(freqs);
bar=waitbar(0, sprintf('Progress: 0/%d', max_i));

for wf_i=1:length(waveforms)
    for dtm_i=1:length(dtheta_maxes)
        for freq_i=1:length(freqs)
            progress=i/max_i;
            waitbar(progress, bar, sprintf('Progress: %d/%d', i, max_i))
            
            freq=freqs(freq_i);
            dtm=dtheta_maxes(dtm_i);
            wf=waveforms{wf_i};
            
            parfor run=1:runs
                [d1, tvec] = test_bench(40000, dtm, freq, wf);

                rhos(wf_i,dtm_i,freq_i,run)=sum(d1.model.rho(d1.nstart:end)) / (d1.model.nsamples - d1.nstart + 1);
                energies(wf_i,dtm_i,freq_i,run)=sum(d1.trg)/d1.model.fs;
            end

            i=i+1;
        end
    end
end

toc

close(bar)

%% Plots

rho_heatmaps = cell(size(waveforms));
energy_heatmaps = cell(size(waveforms));

for i=1:length(waveforms)
    wf_rhos=mean(squeeze(rhos(i,:,:,:)), 3);
    wf_energies=mean(squeeze(energies(i,:,:,:)), 3);

    figure
    px(1)=subplot(2,1,1);
    rho_heatmaps{i} = heatmap(freqs, dtheta_maxes, wf_rhos);
    ylabel('Maximum Waveform Pertubation');
    xlabel('Frequency');
    title(sprintf('%s Wave effect on Synchonization', waveform_names{i}));
    
    px(2)=subplot(2,1,2 );
    energy_heatmaps{i} = heatmap(freqs, dtheta_maxes, wf_energies);
    ylabel('Maximum Waveform Pertubation');
    xlabel('Frequency');
    title(sprintf('Energy Usage for %s Waves', waveform_names{i}));
end

all_rhos = reshape(rhos, [], 1);
all_energies = reshape(energies, [], 1);

for i=1:length(waveforms)
    rho_heatmaps{i}.ColorLimits = [0 max(all_rhos)];
    energy_heatmaps{i}.ColorLimits = [0 max(all_energies)];
end

% Rhos vs Energies Scatter
figure
scatter(all_rhos, all_energies);
xscale('log')
yscale('log')
ylabel('Synchronization');
xlabel('Energy Usage');


%% Save Data
filename = sprintf('out/rho_%s.csv', datestr(now,'mm-dd-yyyy HH-MM'));
writematrix(rhos, filename)
filename = sprintf('out/energy_%s.csv', datestr(now,'mm-dd-yyyy HH-MM'));
writematrix(energies, filename)