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
            %
            % Returns
            % -------
            % - kernel: (neuron x feature x time x latency) double array
            
            folder = this.path.datastore;
            
            kernel = randn(41, 25, 1081, 27);
        end
        
        function rf = getRf(this)
            % Get `receptive field` center for each neuron in dva
            %
            % Returns
            % -------
            % - rf: (neuron x 2) double array
            
            folder = this.path.datastore;
            kernel = this.data.kernel;
            
            rf = randn(41, 2);
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
            
            kernel = this.data.kernel;
            nna = randn(41, 3);
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

