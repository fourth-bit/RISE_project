classdef dbs <handle
   
    properties
        
        %Time index to start stimulation
        nstart double;
        %Maximum frequency of stimulation
        fm double;
        %Charge time to trigger
        ch_tt double;
        %Maximum charge on plate
        Qmax double;
        
        %Number of recording electrodes
        nelec_rec double;
        %Number of stimulating electrodes
        nelec_stim double;
        
        %Charges on plates
        Q double;
        %Triggers from stimulation
        trg double;
        %Triggers for each contact
        trg_l double;
        %Stim artefact from each contact
        stim_l double;
        
        %Channel oscillations (i.e. measurement from electrodes)
        channel double;
        %Model data
        model;
        %Stimulation parameters
        stim_par struct;
        %Matrix to obtain effect of all populations at an electrode.
        D double;
        %Matrix to obtain effect of all electrodes at a population.
        TD double;
        %Stimulation strategy.
        sfunc;
        %Stimulation artefact function.
        sa_func;

        % Waveform function
        wf_func;
        
    end
    
    methods
        
        
        function plot_trg_l(obj)
           
            obj.model.disp_splash();
            px=zeros(1,obj.nelec_stim);
            
            for j=1:obj.nelec_stim
                
                px(j)=subplot(obj.nelec_stim,1,j );
                plot(obj.model.tvec,obj.trg_l(j,:));
                xlabel('Time (s)');
                ylabel(['Trigger ',num2str(j)]);
                
            end
            
            linkaxes(px,'x');
            
        end
        
        
        function plot_channels(obj)
           
            obj.model.disp_splash();
            px=zeros(1,obj.nelec_rec);
            
            for j=1:obj.nelec_rec
                
                px(j)=subplot(obj.nelec_rec,1,j );
                plot(obj.model.tvec,obj.channel(j,:));
                xlabel('Time (s)');
                ylabel(['Channel ',num2str(j)]);
                
            end
            
            linkaxes(px,'x');
            
        end
        
        function plot_osc_trg(obj)
           
            obj.model.disp_splash();
            
            px=zeros(1,2);
            
            px(1)=subplot(2,1,1 );
            plot(obj.model.tvec,obj.model.osc);
            %xlabel('Time');
            %ylabel('Osc');
            ylabel('F(t)');
            
            px(2)=subplot(2,1,2 );
            plot(obj.model.tvec,obj.trg);
            xlabel('Time (s)');
            %ylabel('Trigger');
            ylabel('Stimulation');
            
            linkaxes(px,'x');
            
            
        end

        function plot_osc_rho(obj)
           
            obj.model.disp_splash();
            
            px=zeros(1,2);
            
            px(1)=subplot(2,1,1 );
            plot(obj.model.tvec,obj.model.osc);
            %xlabel('Time');
            %ylabel('Osc');
            ylabel('F(t)');
            
            px(2)=subplot(2,1,2 );
            plot(obj.model.tvec,smooth(obj.model.rho,300));
            xlabel('Time (s)');
            %ylabel('Trigger');
            ylabel('$\rho$', 'interpreter', 'latex');
            
            linkaxes(px,'x');
            
            
        end
        
        function plot_osc_energy(obj)
           
            obj.model.disp_splash();
            
            px=zeros(1,2);
            
            px(1)=subplot(2,1,1 );
            plot(obj.model.tvec,obj.model.osc);
            xlabel('Time (s)');
            ylabel('Osc');
            
            px(2)=subplot(2,1,2 );
            plot(obj.model.tvec,cumsum(obj.trg)/obj.model.fs);
            xlabel('Time (s)');
            ylabel('Energy');
            
            linkaxes(px,'x');
            
        end
        
        function plot_rho_energy(obj)
           
            obj.model.disp_splash();
            
            px=zeros(1,2);
            px(1)=subplot(2,1,1 );
            plot(obj.model.tvec,smooth(obj.model.rho,300));
            xlabel('Time (s)');
            ylabel('$\rho$','interpreter','latex');
            
            px(2)=subplot(2,1,2 );
            plot(obj.model.tvec,cumsum(obj.trg)/obj.model.fs);
            xlabel('Time (s)');
            ylabel('Energy');
            
            linkaxes(px,'x');
            
        end
        
        function update(obj)
           
            % Run Kuramoto
            obj.model.update();
            
            % Grab Kuramoto outputs
            X=zeros(obj.model.npop+obj.nelec_stim,1);
            X(1:obj.model.npop)=obj.model.osc_s(:,obj.model.n);
            % Same action of setting the channels of stimulate
            obj.channel(:,obj.model.n)=(obj.D)*X;
            
            
        end
        
        
        function Q=stimulate(obj)
            
            n=obj.model.n;
            
            if obj.model.n_t<obj.nstart
                Q=zeros(1,obj.nelec_stim);
                obj.Q=Q;
                return;
            end
            
            % Sampling rate / maximum freqency: Tells us how many
            % iterations until we reach charge time to trigger
            nm=round(obj.model.fs/obj.fm);

            % Stimulation function enters here
            stimQ=obj.sfunc(obj);
            
            % Create the waveform from the function
            obj.Q=stimQ * diag(obj.wf_func(obj.ch_tt/nm));
            %obj.Q = stimQ;

            % For anything that hasn't reached charge time to trigger, set
            % the current to 0
            stimQ(obj.ch_tt<nm)=0;
            
            % Show the current triggers
            obj.trg_l(:,n)=obj.Q';
            
            % If it still triggers, set the charge time to zero again
            obj.ch_tt( stimQ>0 )=1;
            % Let everything set to 0 charge up
            obj.ch_tt( stimQ==0 ) =obj.ch_tt( stimQ==0 )+1;
            % Limit things that have reached the threshold trigger
            obj.ch_tt( obj.ch_tt>nm )=nm;
            
            % Use the fired electrodes to change the population Voltage
            obj.model.Vpp=obj.TD*obj.Q';
            % Sum up total triggers
            obj.trg(n)=sum(abs(obj.Q));
            
            % This is the stimulation artefact
            obj.stim_l(:,n)=obj.sa_func(obj);
            
            % Run the kuramoto ODEs, collect outputs in X
            X=zeros(obj.model.npop+obj.nelec_stim,1);
            X(1:obj.model.npop)=obj.model.osc_s(:,n);
            
            % Add in the artefacts
            X(obj.model.npop+1:end)=obj.stim_l(1:obj.nelec_stim,n);
            % Set the channel recordings
            obj.channel(:,n)=(obj.D)*X;
            
        end
        
        function simulate(obj)
           
            obj.model.disp_splash();
            
            initj=obj.model.n+1;
            maxj=obj.model.nsamples;
            bar=waitbar(0, 'Progress: 0%');

            for j=initj:maxj
                progress=(j - initj)/(maxj-initj);
                waitbar(progress, bar, sprintf('Progress: %d%%', floor(progress*100)))

                obj.update();
                obj.stimulate();
                
            end
           
            close(bar);
        end
        
        function initialise_parameters(obj)
        
            obj.TD=zeros(obj.model.npop,obj.nelec_stim);
            obj.D=zeros(obj.nelec_rec,obj.nelec_stim+obj.model.npop);
            
        end
            
        function initialise_arrays(obj)
           
            obj.ch_tt=zeros(1,obj.nelec_stim)+round(obj.model.fs/obj.fm);
            obj.Q=zeros(1,obj.nelec_stim);
            obj.channel=zeros(obj.nelec_rec,obj.model.nsamples);
            obj.trg=zeros(1,obj.model.nsamples);
            obj.trg_l=zeros(obj.nelec_stim,obj.model.nsamples);
            
        end
        
        function obj = dbs(nelec_rec,nelec_stim,fm,qmax,model,sfunc)
            
            obj.model=model;
            obj.nelec_rec=nelec_rec;
            obj.nelec_stim=nelec_stim;
            obj.fm=fm;
            obj.sfunc=sfunc;
            obj.Qmax=qmax;
            obj.stim_par=struct();
            obj.initialise_arrays();
            obj.initialise_parameters();
            
        end
        
       
    end
    
end

