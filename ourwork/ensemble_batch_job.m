%Ensembling Everything into One Script
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

function ensemble_batch_job(nslots)
    %% Setup
    % redirects ~/.matlab PCT temp files to TMPDIR on the compute
    % node to avoid inter-node (compute node <--> login node) I/O
    
    myCluster = parcluster('local'); % cores on compute node to be "local"
    if getenv('ENVIRONMENT')    % true if this is a batch job
      myCluster.JobStorageLocation = getenv('TMPDIR');  % points to TMPDIR
    end
    parpool(myCluster, nslots)    % for MATLAB R2014a or newer

    addpath(genpath('..'));
    
    %waveform_names={'Square', 'Expontial Rise', 'Exponential Decay', 'Sinusoidal', 'Triangular'};
    %waveforms={@wf_none, @wf_exp, @wf_inv_exp, @wf_sinusoid, @wf_triangle};
    waveform_names={'Square' 'Sinusoidal'};
    waveforms={@wf_none, @wf_sinusoid};
    
    dtheta_maxes=linspace(0, 0.0005, 5);
    freqs=linspace(90, 150, 4);
    
    runs=nslots;
    rhos=zeros(length(waveforms), length(dtheta_maxes), length(freqs), runs);
    energies=zeros(length(waveforms), length(dtheta_maxes), length(freqs), runs);
    
    
    %% Run Ensemble
    tic
    
    % For reproduction
    rng('default')
    
    i=0;

    for wf_i=1:length(waveforms)
        for dtm_i=1:length(dtheta_maxes)
            for freq_i=1:length(freqs)                
                freq=freqs(freq_i);
                dtm=dtheta_maxes(dtm_i)*2*pi;
                wf=waveforms{wf_i};
                
                tic
                parfor run=1:runs
                    rng(run);

                    [d1, tvec] = test_bench(40000, dtm, freq, wf);
    
                    rhos(wf_i,dtm_i,freq_i,run)=sum(d1.model.rho(d1.nstart:end)) / (d1.model.nsamples - d1.nstart + 1);
                    energies(wf_i,dtm_i,freq_i,run)=sum(d1.trg)/d1.model.fs;
                end
                toc

                filename=sprintf('group7_out_rho-cp%d_%s.csv', i, datestr(now,'mm-dd-yyyy HH-MM'));
                writematrix(rhos, filename);
                filename=sprintf('group7_out_energy-cp%d_%s.csv', i, datestr(now,'mm-dd-yyyy HH-MM'));
                writematrix(energies, filename);
    
                i=i+1;
            end
        end
    end
    
    toc
    
    %% Save Data
    filename = sprintf('group7_out_rho_%s.csv', datestr(now,'mm-dd-yyyy HH-MM'));
    writematrix(rhos, filename)
    filename = sprintf('group7_out_energy_%s.csv', datestr(now,'mm-dd-yyyy HH-MM'));
    writematrix(energies, filename)

    delete(gcp);
end