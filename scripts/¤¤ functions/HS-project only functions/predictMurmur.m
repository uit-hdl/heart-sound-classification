function Y = predictMurmur(id,net,varargin)
% takes a a set of four recordings (order: aortic, pulmonic, tricuspic,
% mitral) and a trained neural network and outputs an array Y of murmur
% grade prediction.

%% example input
% load 'networksTrainingSet_valStop.mat'
% id = HSdata.id(3);
% % ground truth: [2, 0.5, 0, 0]
%% optional argumets
noise_index  = [0,0,0,0];
N_downSample       = 20;
N_cyclesPerSegment = 4;
N_cycleOverlap     = 2;
N_segmentsPerHSrec = 6;
N_downSample_ACF   = 35;
MFCC_sz = [13,200];
HMMpar = [];

p = inputParser;
addOptional(p, 'noise_index', noise_index, @(x) isnumeric(x));
addOptional(p, 'N_downSample', N_downSample, @(x) isnumeric(x));
addOptional(p, 'N_cyclesPerSegment', N_cyclesPerSegment, @(x) isnumeric(x));
addOptional(p, 'N_cycleOverlap', N_cycleOverlap, @(x) isnumeric(x));
addOptional(p, 'N_segmentsPerHSrec', N_segmentsPerHSrec, @(x) isnumeric(x));
addOptional(p, 'N_downSample_ACF', N_downSample_ACF, @(x) isnumeric(x));
addOptional(p, 'HMMpar', HMMpar, @(x) isnumeric(x));
parse(p,varargin{:});

noise_index  = p.Results.noise_index;
N_downSample       = p.Results.N_downSample;
N_cyclesPerSegment = p.Results.N_cyclesPerSegment;
N_cycleOverlap     = p.Results.N_cycleOverlap;
N_segmentsPerHSrec = p.Results.N_segmentsPerHSrec;
N_downSample_ACF   = p.Results.N_downSample_ACF;
HMMpar = p.Results.HMMpar;

%% body

if isempty(HMMpar)
    load('HMMpar.mat','HMMpar')
end

% *** segmentation ***
segmentation = ngbrSegment(id, N_downSample, N_downSample_ACF, HMMpar);

% container for murmur predictions:
Y = zeros(4,1);

for aa=1:4
    % container for activations for position aa:
    A = zeros(N_segmentsPerHSrec,1);
    
    if noise_index(aa)==0
        
        % *** get and preproces audio ***
        [x0,fs0] = wav2TS(id,aa);
        x  = schmidt_spike_removal(x0, fs0);
        x  = downsample(x,N_downSample);
        fs = floor(fs0/N_downSample);
        
        % compute lines where diastole ends (segmentation lines):
        [~,~,segLines] = getCompactReprOfStates(segmentation{aa});
        
        % in case there are fewer cycles than the number requested per segment:
        N_cyclesAvailable = numel(segLines) - 1;
        N_cyclesPerSegment_actually  = min(N_cyclesPerSegment, N_cyclesAvailable);
        
        % get beginning and end of each segment:
        [L,~] = getSegments(segLines, N_cycleOverlap, N_cyclesPerSegment_actually, N_segmentsPerHSrec);

        % *** cycle through the segments and compute MFCC input ***
        for k=1:N_segmentsPerHSrec
            % get the segment:
            xk = x(L(k,1):L(k,2));
            % get compact representation of signal in time-frequency domain:
            Mk = getMFCC(xk,fs);
            % normalize MFCC:
            Mk = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
            % resize:
            Mk = imresize(Mk, MFCC_sz);
            % predict and store prediction:
            A(k) = predict(net, Mk);
        end
        % take median to get predictions
        Y(aa) = median(A);
        
    end
end

end