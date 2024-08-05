%Plotting everything from the ensemble

waveform_names={'Rectangular', 'Sinusoidal', 'Exponential Decay'};
waveforms={@wf_none, @wf_sinusoid, @wf_inv_exp};

dtheta_maxes=linspace(0, 0.0005, 5);
freqs=linspace(90, 150, 4);

runs=28;

rhos=readmatrix('out/complete_rho_08-01.csv');
energies=readmatrix('out/complete_energy_08-01.csv');

rhos = reshape(rhos, [length(waveforms), length(dtheta_maxes), length(freqs), runs]);
energies = reshape(energies, [length(waveforms), length(dtheta_maxes), length(freqs), runs]);

rho_heatmaps = cell(size(waveforms));
energy_heatmaps = cell(size(waveforms));

max_rho=0;
min_rho=0;
max_energy=0;

for i=1:length(waveforms)
    wf_rhos=mean(squeeze(rhos(i,:,:,:)), 3);
    wf_energies=mean(squeeze(energies(i,:,:,:)), 3);

    max_rho=max(max_rho, max(wf_rhos, [], 'all'));
    min_rho=min(min_rho, min(wf_rhos, [], 'all'));
    max_energy=max(max_energy, max(wf_energies, [], 'all'));

    figure
    px(1)=subplot(2,1,1);
    rho_heatmaps{i} = heatmap(freqs, dtheta_maxes, wf_rhos);
    colormap(flipud(summer));
    ylabel('Maximum Waveform Pertubation');
    xlabel('Frequency');
    title(sprintf('%s Wave Effect on Synchonization', waveform_names{i}));
    fontsize(20, 'points');
    fontname('Calibri');  
    
    px(2)=subplot(2,1,2);
    energy_heatmaps{i} = heatmap(freqs, dtheta_maxes, wf_energies);
    colormap(energy_heatmaps{i}, flipud(winter));
    ylabel('Maximum Waveform Pertubation');
    xlabel('Frequency');
    title(sprintf('Energy Usage for %s Waves', waveform_names{i}));
    fontsize(20, 'points');
    fontname('Calibri');
end

for i=1:length(waveforms)
    rho_heatmaps{i}.ColorLimits = [min_rho max_rho];
    energy_heatmaps{i}.ColorLimits = [0 max_energy];
end

%% Rhos vs Energies Scatter
all_rhos = reshape(mean(rhos, 4), [], 1);
all_energies = reshape(mean(energies, 4), [], 1);

figure
[f, goodness]=fit(all_energies, all_rhos, 'exp1');
scatter(all_energies, all_rhos, 100, 'filled');
ylimits = ylim();

hold on

h=plot(f);
set(h, 'LineWidth', 1.5);
ylabel('Synchronization');
xlabel('Energy Usage');
ylim(ylimits);
fontsize(20, 'points');
fontname('Calibri');

fn_str=sprintf('%.2fe^{%.2fx}', f.a, f.b);
legend({'', fn_str});

hold off

%% Some analysis
rho_significance = zeros(length(waveforms), length(waveforms), length(dtheta_maxes), length(freqs));
energy_significance = zeros(length(waveforms), length(waveforms), length(dtheta_maxes), length(freqs));

for i=1:length(waveforms)
    for j=1:length(waveforms)
        for dtm_i=1:length(dtheta_maxes)
            for freq_i=1:length(freqs)
                [~,rho_significance(i,j,dtm_i,freq_i)]=ttest(rhos(i, dtm_i, freq_i,:),rhos(j, dtm_i, freq_i,:));
                [~,energy_significance(i,j,dtm_i,freq_i)]=ttest(energies(i, dtm_i, freq_i,:),energies(j, dtm_i, freq_i,:));
            end
        end
    end
end