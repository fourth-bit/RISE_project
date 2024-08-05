%Simple shape testing setup
%System geometry: spherical
%Frequency distribution: Lorentzian
%Stimulation strategy: ACD

function d1 = run_populations(nelec, npop, eta)
    %% Set model parameters
    %Set number of oscillators per population
    nperpop=600;
    
    %Parameters for natural frequency lorentzian distribution.
    %[width,centre]
    lorentz_param=[0.9615   24.6316];
    lorentzian_centre=lorentz_param(2);
    
    %Simulation time
    T=200;
    %Sampling rate
    fs=400;
    %Noise determined as a multiple of the centre frequency
    sigma_mult=0.05;
    %Off diagonal components of the coupling matrix.
    k_off_diag=0;
    %Diagonal components of the coupling matrix.
    k_on_diag=55;
    
    %Number of samples in simulation.
    nsamples=T*fs;
    %Time domain vector.
    tvec=linspace(0,T,T*fs);
    
    %Generate instantaneous phase response curve
    %zA cosine Fourier coefficients [a0,a1]
    %zB sine Fourier coefficients [b0,b1]
    zA=[4,0];
    zB=[0,-1];
    
    %Integration step
    dt=1/fs;
    
    %Create multi-population Kuramoto configurator structure
    mp_km_model=create_mp_km_model();
    
    %Configure the model
    mp_km_model.T=T;
    mp_km_model.fs=fs;
    mp_km_model.npop=npop;
    mp_km_model.nperpop=nperpop;
    mp_km_model.sigma=(lorentzian_centre)*sigma_mult;
    mp_km_model.kmat=get_kmat(k_on_diag,k_off_diag,mp_km_model.npop);
    mp_km_model.gamma=lorentz_param(1);
    mp_km_model.w0=lorentz_param(2);
    mp_km_model.nat_freq_dist='cauchy';
    mp_km_model.zA=zA;
    mp_km_model.zB=zB;
    
    %Create instance of the multi-population Kuramoto object
    mpk=create_mp_km_obj(mp_km_model);
    
    %% Set DBS parameters
    
    %Generate random system in a unit sphere.
    %Electrodes can stimulate and record
    %P: positions of electrodes.
    %Ppp: positions of populations.
    
    [P,Ppp]=generate_random_system_spherical_fl_ratio(npop,nelec,eta);
    
    %Obtain transformation matrices. D converts vector of population activities
    %into electrode measurements. TD converts vector of charges due to
    %stimulation into a 'stimulation intensity' at a population.
    [D,TD] = get_Dmat_coulombic_single(P,Ppp);
    
    %Set maximum stimulation current (qmax) such that the maximimum perturbation to a
    %particular oscillator in the system in dtheta_max.
    dtheta_max=0.0005*2*pi;
    %Max of Z
    Zmax=1;
    samples_per_fire=round(fs/130) + 1;
    qmax=get_qmax(dtheta_max,Zmax,dt,TD) / sum(wf_none(linspace(0, 1, samples_per_fire)));
    % qmax = qmax / nelec * 3;
    
    %Create dbs model configurator structure. Input name of stimulation
    %strategy function, found in lib/stimulation.
    dbs_model=create_dbs_model('acd');
    %Stimulation start time.
    dbs_model.tstart=15;
    %Maximum stimulation frequency.
    dbs_model.fm=130;
    
    %Configure the model
    dbs_model.model=mpk;
    dbs_model.Qmax=qmax;
    dbs_model.nelec_rec=nelec;
    dbs_model.nelec_stim=nelec;
    dbs_model.D=D;
    dbs_model.TD=TD;
    
    %Set stimulation artefact function, found in model/lib/artefact.
    dbs_model.sa_func=@sa_none;
    %Set amplitude of artefact effect
    dbs_model.sa_amp=0;
    
    dbs_model.wf_func=@wf_none;
    
    %Create configured DBS object
    d1=create_dbs_obj(dbs_model);
    
    %% Run simulation
    
    d1.simulate_no_bar();
end

%% Simulation
% We have a 3d parameter space: nelec, npop, and eta
etas=[0.1 0.3 0.8];
% Hold nelec at 4
nelec=4;
npops=[4 6 8];

runs=32;
pop2_dbs=cell(length(etas), length(npops), runs);
rho_avg_start=20;
nrho_avg_start=rho_avg_start*400+1;

tic
for eta_i=1:length(etas)
    eta=etas(eta_i);
    for npop_i=1:length(npops)
        npop=npops(npop_i);

        tic
        parfor run=1:runs
            pop2_dbs{eta_i, npop_i, run} = run_populations(nelec, npop, eta);
        end
        toc

        fprintf('Completed eta=%g, npop=%d\n', eta, npop);
    end
end
toc


%% Process Data

% Threshold rho? I.e. a metric like samples above rho? 
% Standard deviation?

raw_data=zeros(length(etas), length(npops), runs);
avgs=zeros(length(etas), length(npops));
stdevs=zeros(length(etas), length(npops));

for eta_i=1:length(etas)
    integrals=zeros(length(npops), runs);

    for npop_i=1:length(npops)
        for run_i=1:runs
            d1=pop2_dbs{eta_i, npop_i, run_i};
            d1_rho=d1.model.rho(nrho_avg_start:end);
            integral = sum(d1_rho) / 400;
            integrals(npop_i, run_i)=integral;
        end
    end

    avgs(eta_i,:) = mean(integrals, 2);
    stdevs(eta_i,:) = std(integrals, 0, 2);
    raw_data(eta_i,:,:) = integrals;
end

%% Plots
eta1=squeeze(raw_data(1,:,:));
eta2=squeeze(raw_data(2,:,:));
eta3=squeeze(raw_data(3,:,:));

figure
subplot(2, 1, 1)
p = boxplot(eta1.', npops);
ylabel('Total Syncrhony')
% Already ran the t-tests (ttest2)
sigstar({[1,2], [2,3], [1,3]}, [0.00065 0.0065 0.00065])
fontsize(25, 'points')
set(p, 'LineWidth', 1.5)
title('\eta = 0.1')

subplot(2, 1, 2)
p = boxplot(eta3.', npops);
xlabel('Number of Populations')
ylabel('Total Synchrony')
% Already ran the t-tests (ttest2)
sigstar({[1,2], [1,3]}, [0.0015 0.0184])
fontsize(25, 'points')
set(p, 'LineWidth', 1.5)
title('\eta = 0.8')

%% Save Data
filename=sprintf('out/changed_population_size-%s.csv', datestr(now,'mm-dd-yyyy HH-MM'));
writematrix(raw_data, filename);