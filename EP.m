classdef EP < handle
    % EffectPropagation: Propagation of perisaccadic effect across V4/MT
    % population
    
    properties
        % - path
        %   - assets: assets directory
        %   - datastore: recorded and trained data directory
        %   - data: filename of input data
        %   - x: filename of `x` features
        %   - y: filename of `y` features
        %   - rho
        %       - d1Raw: filename of correlation between `d1` and `raw`
        %       - d1Nna: filename of correlation between `d1` and `nna`
        % - population: number of neurons in population
        % - data
        %   - kernel: (neuron x feature x time x latency) double
        %   array
        %   - rf: (neuron x 2) double array. (x, y) `receptive field`
        %   center for each neuron at fixation time in dva
        %   - fp1: (neuron x 2) double array. (x, y) `fixation point 1` for
        %   each neuron in dva
        %   - fp2: (neuron x 2) double array. (x, y) `fixation point 2` for
        %   each neuron in dva
        % - idx
        %   - rf
        %       - a: amplitude
        %       - mx: mu_x
        %       - my: mu_y
        %       - sx: sigma_x
        %       - sy: sigma_y
        %       - r: rho
        %       - gx: gamma_x
        %       - gy: gamma_y
        %   - ff
        %       - a
        %       - mx
        %       - my
        %       - sx
        %       - sy
        %       - r
        %       - gx
        %       - gy
        %   - st
        %       - a
        %       - mx
        %       - my
        %       - sx
        %       - sy
        %       - r
        %       - gx
        %       - gy
        %   - time: from saccadic onset (ms)
        %       - fix:      -inf:-100
        %       - pre:      -100:0
        %       - sacon:    0
        %       - sacoff:   30
        %       - post:     0:100
        %       - peri:     -100:100
        %       - st:       100:inf
        % - x
        %   - d1: (neuron x 1) double vector. distance between `RF` and
        %   `FP1` for each neuron
        %   - d2: (neuron x 1) double vector. distance between `RF` and
        %   `FP2` for each neuron
        %   - sd1: (neuron x 1) double vector. signed distance (twoard
        %   `FP2` is positive) between `RF` and `FP1` for each neuron
        %   - sd2: (neuron x 1) double vector. signed distance (twoard
        %   `FP2` is positive) between `RF` and `FP2` for each neuron
        % - y
        %   - raw: (neuron x feature x time x latency) double array same as
        %   kernel of each neuron
        %   - nna: (neuron x 3) double matrix. extract `nna1`, `nna2` and
        %   `nna3` features for each neuron
        path
        population
        data
        idx
        x
        y
    end
    
    % Constructor
    methods
        function this = EP()
            % Constructor
            
            % path
            this.initPath();
            
            % population
            this.population = 41;
            
            % data
            disp('Data');
            tic();
            this.initData();
            toc();
            
            % x
            disp('X features');
            tic();
            this.initX();
            toc();
            
            % y
            disp('Y features');
            tic();
            this.initY();
            toc();
        end
    end
    
    % Path
    methods
        function initPath(this)
            % assets
            assets = './assets';
            
            this.path.assets = assets;
            
            if ~exist('./assets', 'dir')
                mkdir(assets);
            end
            
            % datastore
            % <kaiser-amir>
            %             this.path.datastore = 'path/to/datastore/';
            this.path.datastore = 'K:\Barfak\scdata';

            % </kaiser-amir>
            
            % data
            this.path.data = fullfile(assets, 'data.mat');
            
            % x
            this.path.x = fullfile(assets, 'x.mat');
            
            % y
            this.path.y = fullfile(assets, 'y.mat');
            
            % rho
            % - d1, raw
            this.path.rho.d1Raw = fullfile(assets, 'd1_raw.mat');
            % - d1, nna
            this.path.rho.d1Nna = fullfile(assets, 'd1_nna.mat');
        end
    end
    
    % Data
    methods
        function initData(this)
            % `data` file contains `kernel`, `rf`, `fp1` and `fp2`
            filename = this.path.data;
            
            % load saved `data` file
            if exist(filename, 'file')
                load(filename, 'kernel', 'rf', 'fp1', 'fp2');
            end
            
            % `kernel`
            if ~exist('kernel', 'var')
                kernel = this.getKernel();
                save(filename, 'kernel');
            end
            this.data.kernel = kernel;
            
            % `rf`
            if ~exist('rf', 'var')
                rf = this.getRf();
                save(filename, 'rf', '-append');
            end
            this.data.rf = rf;
            
            % `fp1`
            if ~exist('fp1', 'var')
                fp1 = this.getFp1();
                save(filename, 'fp1', '-append');
            end
            this.data.fp1 = fp1;
            
            % `fp2`
            if ~exist('fp2', 'var')
                fp2 = this.getFp2();
                save(filename, 'fp2', '-append');
            end
            this.data.fp2 = fp2;
        end
        
        % <kaiser>
        function kernel = getKernel(this)
            % Get `kernel` for each neuron
            % Returns
            % - kernel: (neuron x feature x time x latency) double array
            kernel = nan(this.population,25,1081,27);
            dblst = load('J:\MATLAB\MT Gaussian Modelling\Paper Figures\db\make_the_main_lists.mat');
            for int_neuron = 1:size(dblst.probe_list,1)
                sess = dblst.probe_list(int_neuron,1);
                chnl = dblst.probe_list(int_neuron,2);
                mdl = load(['J:\MATLAB\MT Gaussian Modelling\Effect Studies\modulation problem\Gaussian\models\',num2str(sess),'_',num2str(chnl),'\fold_0.mat']);                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                try
                    scd = load([this.path.datastore,'scdata_',num2str(dblst.probe_list(int_neuron,1)-20000000),'_',strcat(num2str((dblst.probe_list(int_neuron,2)-mod(dblst.probe_list(int_neuron,2),10))/10),num2str(mod(dblst.probe_list(int_neuron,2),10))),'_1.mat'],'tdata');
                    [fp_x,fp_y] = pol2cart(scd.tdata.params.Fix1_theta/180*pi,scd.tdata.params.Fix1_radius);
                    [ct_x,ct_y] = pol2cart(scd.tdata.params.prob_theta/180*pi,scd.tdata.params.probe_radius);
                    probTweak_x = scd.tdata.params.probTweak_x;
                    probTweak_y = scd.tdata.params.probTweak_y;
                    probsize = scd.tdata.params.probsize;
                catch
                    fp_x = nan; fp_y = nan;
                    ct_x = nan; ct_y = nan;
                    probTweak_x = nan; probTweak_y = nan;
                    probsize = nan;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set_of_coefficients = mdl.nProfile.set_of_kernels.gas.cof;
                modified_set_of_coefficients = nan(size(set_of_coefficients));
                % RF part
                modified_set_of_coefficients(:,:,1) = set_of_coefficients(:,:,1);
                modified_set_of_coefficients(:,:,2) = fp_x + ct_x + (set_of_coefficients(:,:,2) - (probsize+1)/2)*probTweak_x; % mu_x
                modified_set_of_coefficients(:,:,3) = fp_y + ct_y + (set_of_coefficients(:,:,3) - (probsize+1)/2)*probTweak_y; % mu_y
                modified_set_of_coefficients(:,:,4) = set_of_coefficients(:,:,4).*probTweak_x; % sigma_x
                modified_set_of_coefficients(:,:,5) = set_of_coefficients(:,:,5).*probTweak_y; % sigma_y
                modified_set_of_coefficients(:,:,6) = set_of_coefficients(:,:,6); % rho
                modified_set_of_coefficients(:,:,7) = set_of_coefficients(:,:,7)./probTweak_x; % gamma_x
                modified_set_of_coefficients(:,:,8) = set_of_coefficients(:,:,8)./probTweak_y; % gamma_y
                % FF part
                modified_set_of_coefficients(:,:,09) = set_of_coefficients(:,:,09);
                modified_set_of_coefficients(:,:,10) = fp_x + ct_x + (set_of_coefficients(:,:,10) - (probsize+1)/2)*probTweak_x; % mu_x
                modified_set_of_coefficients(:,:,11) = fp_y + ct_y + (set_of_coefficients(:,:,11) - (probsize+1)/2)*probTweak_y; % mu_y
                modified_set_of_coefficients(:,:,12) = set_of_coefficients(:,:,12).*probTweak_x; % sigma_x
                modified_set_of_coefficients(:,:,13) = set_of_coefficients(:,:,13).*probTweak_y; % sigma_y
                modified_set_of_coefficients(:,:,14) = set_of_coefficients(:,:,14); % rho
                modified_set_of_coefficients(:,:,15) = set_of_coefficients(:,:,15)./probTweak_x; % gamma_x
                modified_set_of_coefficients(:,:,16) = set_of_coefficients(:,:,16)./probTweak_y; % gamma_y
                % ST part
                modified_set_of_coefficients(:,:,17) = set_of_coefficients(:,:,17);
                modified_set_of_coefficients(:,:,18) = fp_x + ct_x + (set_of_coefficients(:,:,18) - (probsize+1)/2)*probTweak_x; % mu_x
                modified_set_of_coefficients(:,:,19) = fp_y + ct_y + (set_of_coefficients(:,:,19) - (probsize+1)/2)*probTweak_y; % mu_y
                modified_set_of_coefficients(:,:,20) = set_of_coefficients(:,:,20).*probTweak_x; % sigma_x
                modified_set_of_coefficients(:,:,21) = set_of_coefficients(:,:,21).*probTweak_y; % sigma_y
                modified_set_of_coefficients(:,:,22) = set_of_coefficients(:,:,22); % rho
                modified_set_of_coefficients(:,:,23) = set_of_coefficients(:,:,23)./probTweak_x; % gamma_x
                modified_set_of_coefficients(:,:,24) = set_of_coefficients(:,:,24)./probTweak_y; % gamma_y
                % c part
                modified_set_of_coefficients(:,:,25) = set_of_coefficients(:,:,25); % c
                % permute, and insert
                t_modified_set_of_coefficients = permute(modified_set_of_coefficients,[3 1 2]);
                kernel(int_neuron,:,:,:) = t_modified_set_of_coefficients;
                clearvars modified_set_of_coefficients t_modified_set_of_coefficients
                clearvars fp_x fp_y ct_x ct_y probTweak_x probTweak_y probsize mdl sess chnl
            end
        end
        
        function rf = getRf(this)
            % Get `receptive field` center for each neuron in dva
            % Returns
            % - rf: (neuron x 2) double array
            rf = nan(this.population,2);
            dblst = load('J:\MATLAB\MT Gaussian Modelling\Paper Figures\db\make_the_main_lists.mat');
            for int_neuron = 1:size(dblst.probe_list,1)
                rf_in_probe = dblst.probe_list(int_neuron,3:4);
                try
                    scd = load([this.path.datastore,'scdata_',num2str(dblst.probe_list(int_neuron,1)-20000000),'_',strcat(num2str((dblst.probe_list(int_neuron,2)-mod(dblst.probe_list(int_neuron,2),10))/10),num2str(mod(dblst.probe_list(int_neuron,2),10))),'_1.mat'],'tdata');
                    [fp_x,fp_y] = pol2cart(scd.tdata.params.Fix1_theta/180*pi,scd.tdata.params.Fix1_radius);
                    [ct_x,ct_y] = pol2cart(scd.tdata.params.prob_theta/180*pi,scd.tdata.params.probe_radius);
                    rf_in_dva = [ ...
                        fp_x + ct_x + (rf_in_probe(1) - (scd.tdata.params.probsize+1)/2)*scd.tdata.params.probTweak_x , ...
                        fp_y + ct_y + (rf_in_probe(2) - (scd.tdata.params.probsize+1)/2)*scd.tdata.params.probTweak_y ];
                    rf(int_neuron,:) = rf_in_dva;
                catch
                    rf(int_neuron,:) = [nan nan];
                end
            end
        end
        % </kaiser>
        
        % <amir>
        function fp1 = getFp1(this)
            % Get `fixation point 1` for each neuron in dva
            %
            % Returns
            % -------
            % - fp1: (neuron x 2) double array
            
            %             folder = this.path.datastore;
            fp1 = zeros(41,2);
            %             fp1 = randn(41, 2);
        end
        
        function fp2 = getFp2(this)
            % Get `fixation point 2` for each neuron in dva
            %
            % Returns
            % -------
            % - fp2: (neuron x 2) double array
            
            folder = this.path.datastore;
            info = load('J:\Projects\m190124_spatial_sensitivity\m190205_state_based_neurometric_y\make_the_main_lists.mat');
            fp2 = nan(size(info.nuron_list,1),2);
            for nn = 1:size(info.nuron_list,1)
                clear xx
                id = num2str(info.nuron_list(nn,1)-20000000);
                if(info.nuron_list(nn,2)>=10)
                    ch = num2str(info.nuron_list(nn,2));
                else
                    ch = ['0' num2str(info.nuron_list(nn,2))];
                end
                un = '1';
                xx = load([folder '\scdata_' id '_' ch '_' un '.mat'],'trial_info');
                fp2(nn,:) = xx.trial_info.saccade_target;
            end
        end
        % </amir>
    end
    
    % Idx
    methods
        function initIdx(this)
            % Initialize `path`
            
            % rf
            this.idx.rf.a   =    1;
            this.idx.rf.mx  =    2;
            this.idx.rf.my  =    3;
            this.idx.rf.sx  =    4;
            this.idx.rf.sy  =    5;
            this.idx.rf.r   =    6;
            this.idx.rf.gx  =    7;
            this.idx.rf.gy  =    8;
            
            % ff
            this.idx.ff.a   =    9;
            this.idx.ff.mx  =    10;
            this.idx.ff.my  =    11;
            this.idx.ff.sx  =    12;
            this.idx.ff.sy  =    13;
            this.idx.ff.r   =    14;
            this.idx.ff.gx  =    15;
            this.idx.ff.gy  =    16;
            
            % st
            this.idx.st.a   =   17;
            this.idx.st.mx  =   18;
            this.idx.st.my  =   19;
            this.idx.st.sx  =   20;
            this.idx.st.sy  =   21;
            this.idx.st.r   =   22;
            this.idx.st.gx  =   23;
            this.idx.st.gy  =   24;
            
            % time
            this.idx.time.fix = [-500:-100]+541;
            this.idx.time.pre = [];
            this.idx.time.sacon = 0+541;
            this.idx.time.sacoff = this.idx.time.sacon + 30;
            this.idx.time.post = [];
            this.idx.time.peri = [];
            this.idx.time.st = [];
            this.idx.time.pa = [-010:+010]+541;
            this.idx.time.ff = [-050:-000]+541;
            this.idx.time.st = [-050:-000]+541;
        end
        
        function v = getParam(neuron, source, name, time, latency)
            % For a neuron, get parameter values of a specified source at 
            % given times and latencies
            %
            % Parameters
            % ----------
            % - source: ['rf', 'ff', 'st']
            %   Source of effect
            % - name: ['a', 'mx', 'my', 'sx', 'sy', 'r', 'gx', 'gy']
            %   Name of parameter
            % - time: integer | integer vector | ['fix', 'pre', 'sacon',
            % 'sacoff', 'post', 'peri', 'st']
            %   Time from saccadic onset (ms)
            % - latency: integer | integer vector
            %
            % Returns
            % -------
            % - v: double array
            %   Values of specified parameter based on `neuron`, `source`, 
            %   `name`, `time` and `latency`
            
            % feature
            feature = this.idx.(source).(name);
            
            % time
            if ~exist('time', 'var')
                time = 1:size(this.kernel, 3);
            end
            if ischar(time)
                time = this.idx.time.(time);
            end
            
            % latency
            if ~exist('latency', 'var')
                latency = 1:size(this.kernel, 4);
            end
            
            % values
            v = this.kernel(neuron, feature, time, latency);
        end
    end
    
    % X features
    methods
        function initX(this)
            % `x` file contains `d1`, `d2`, `sd1` and `sd2`
            filename = this.path.x;
            
            % load saved `x` file
            if exist(filename, 'file')
                load(filename, 'd1', 'd2', 'sd1', 'sd2');
            end
            
            % `d1`
            if ~exist('d1', 'var')
                d1 = this.getD1();
                save(filename, 'd1');
            end
            this.x.d1 = d1;
            
            % `d2`
            if ~exist('d2', 'var')
                d2 = this.getD2();
                save(filename, 'd2', '-append');
            end
            this.x.d2 = d2;
            
            % `sd1`
            if ~exist('sd1', 'var')
                sd1 = this.getSd1();
                save(filename, 'sd1', '-append');
            end
            this.x.sd1 = sd1;
            
            % `sd2`
            if ~exist('sd2', 'var')
                sd2 = this.getSd2();
                save(filename, 'sd2', '-append');
            end
            this.x.sd2 = sd2;
        end
        
        function d1 = getD1(this)
            % Get distance between `RF` and `FP1` for each neuron
            %
            % Returns
            % -------
            % - d1: (neuron x 1) double vector
            
            n = this.population;
            rf = this.data.rf;
            fp1 = this.data.fp1;
            
            d1 = zeros(n, 1);
            for i = 1:n
                d1(i) = EP.dist(fp1(i, :), rf(i, :));
            end
        end
        
        function d2 = getD2(this)
            % Get distance between `RF` and `FP2` for each neuron
            %
            % Returns
            % -------
            % - d2: (neuron x 1) double vector
            
            n = this.population;
            rf = this.data.rf;
            fp2 = this.data.fp2;
            
            d2 = zeros(n, 1);
            for i = 1:n
                d2(i) = EP.dist(fp2(i, :), rf(i, :));
            end
        end
        
        function sd1 = getSd1(this)
            % Get signed distance between `RF` and `FP1` for each neuron
            %
            % Returns
            % -------
            % - sd1: (neuron x 1) double vector
            
            fp1 = this.data.fp1;
            fp2 = this.data.fp2;
            rf = this.data.rf;
            
            d1 = this.x.d1;
            sd1 = zeros(size(d1));
            
            for i = 1:size(d1, 1)
                u = fp2(i, :) - fp1(i, :);
                v = rf(i, :) - fp1(i, :);
                
                if dot(u, v) >= 0
                    sd1 = d1(i);
                else
                    sd1 = -d1(i);
                end
            end
        end
        
        function sd2 = getSd2(this)
            % Get signed distance between `RF` and `FP2` for each neuron
            %
            % Returns
            % -------
            % - sd2: (neuron x 1) double vector
            
            sd1 = this.x.sd1;
            d2 = this.x.d2;
            
            sd2 = d2 * sign(sd1);
        end
    end
    
    % Y features
    methods
        function initY(this)
            % `y` file contains `raw` and `nna`
            filename = this.path.y;
            
            % load saved `y` file
            if exist(filename, 'file')
                load(filename, 'raw', 'nna');
            end
            
            % `raw`
            if ~exist('raw', 'var')
                raw = this.getRaw();
                save(filename, 'raw');
            end
            this.y.raw = raw;
            
            % `nna`
            if ~exist('nna', 'var')
                nna = this.getNna();
                save(filename, 'nna', '-append');
            end
            this.y.nna = nna;
        end
        
        function raw = getRaw(this)
            % Get `raw` features (same as sensitivity values) for each
            % neuron
            %
            % Returns
            % -------
            % - raw: (neuron x feature x time x latency) double array
            
            raw = this.data.kernel;
        end
        
        % <neda>
        function nna = getNna(this)
            % Get `nna` (neda nategh amplitude) features for each neuron
            %
            % Returns
            % -------
            % - nna: (neuron x 3) double matrix
            
            neurons = 1:this.population;
            sm = @EP.summax;
            
            % nna_rf
            % - (a_rf_fix - a_rf_peri) / a_rf_fix
            a_rf_fix = sm(this.getParam(neurons, 'rf', 'a', 'fix'));
            a_rf_peri = sm(this.getParam(neurons, 'rf', 'a', 'peri'));
            nna_rf = (a_rf_fix - a_rf_peri) ./ a_rf_fix;
            
            % nna_ff
            % - (a_ff_fix - a_ff_peri) / a_ff_fix
            a_ff_fix = sm(this.getParam(neurons, 'ff', 'a', 'fix'));
            a_ff_peri = sm(this.getParam(neurons, 'ff', 'a', 'peri'));
            nna_ff = (a_ff_fix - a_ff_peri) ./ a_ff_fix;
            
            % nna_st
            % - (a_rf_fix - a_rf_peri) / a_rf_fix
            a_st_fix = sm(this.getParam(neurons, 'st', 'a', 'fix'));
            a_st_peri = sm(this.getParam(neurons, 'st', 'a', 'peri'));
            nna_st = (a_st_fix - a_st_peri) ./ a_st_fix;
            
            nna = [nna_rf, nna_ff, nna_st];
        end
        % </neda>
    end
    
    % Rho
    methods
        function corrD1Raw(this)
            % Correlation between `d1` and `raw`
            
            % x-feature
            d1 = this.x.d1;
            
            % y-feature
            raw = this.y.raw;
            
            % corr(x, y)
            [s1, s2, s3, s4] = size(raw);
            rho = zeros(s1, s2, s3, s4);
            for i1 = 1:s1
                for i2 = 1:s2
                    for i3 = 1:s3
                        for i4 = 1:s4
                            rho(i1, i2, i3, i4) = corr(d1, squeeze(raw(:, i2, i3, i4)));
                        end
                    end
                end
            end
            
            % save
            filename = this.path.rho.d1Raw;
            save(filename, 'd1', 'raw', 'rho');
        end
        function corrD1Nna(this)
            % Correlation between `d1` and `nna`
            
            % x-feature
            d1 = this.x.d1;
            
            % y-feature
            nna = this.y.nna;
            
            % corr(x, y)
            [s1, s2] = size(nna);
            rho = zeros(s1, s2);
            for i1 = 1:s1
                for i2 = 1:s2
                    rho(i1, i2) = corr(d1, squeeze(nna(:, i2)));
                end
            end
            
            % save
            filename = this.path.rho.d1Raw;
            save(filename, 'd1', 'nna', 'rho');
        end
    end
    
    % Libs
    methods (Static)
        function d = dist(p1, p2)
            % Euclidean distance between two points
            
            d = norm(p1 - p2);
        end
        function v = summax(v)
            % Sum on `time`, max on `latency`
            %
            % Parameters
            % ----------
            % - v: [neuron x feature x time x latency] double array
            %   Input values
            %
            % Returns
            % -------
            % - v: [neuron x 1] double vector
            
            v = squeeze(sum(max(v, [], 4), 3));
        end
    end
    
    % Main
    methods (Static)
        function main()
            % Main
            
            % clear
            close('all');
            clear();
            clc();
            
            % compute correlation
            ep = EP();
            
            % corr(d1, raw)
            disp('corr(d1, raw)');
            tic();
            ep.corrD1Raw();
            toc();
            
            % corr(d1, nna)
            disp('corr(d1, nna)');
            tic();
            ep.corrD1Nna();
            toc();
        end
    end
end