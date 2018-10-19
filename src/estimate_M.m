classdef estimate_M
    %estimate_M Different ways of estimating the prior M(b)
    %   Wrapper with several functions
    
    methods(Static)
        %% Probabilistic with mismatched bBias
        function Y = pHT(w, X)
            %pHT Summary of this function goes here
            %   Detailed explanation goes here
            
            tmp = [{X(:,1); X(:,2); X(:,3)}; num2cell(w')];
            Y = calculate(tmp{:});
            
            function [ pbEv ] = calculate( mEv, N, block, varargin )
                % Set defaults
                optargs = {5 9};
                % Overwrite defaults if provided
                numvarargs = length(varargin);
                optargs(1:numvarargs) = varargin;
                % Place optional args in memorable variable names
                [asym, base] = optargs{:};
                asym = asym + base;   % bBias1 formulated as difference to bBias2 to enforce constraint
                
                % Derive unperturbed momentary block evidence
                [ ~, ~, pbEv, ~ ] = opt_inf.all( mEv.*N, N, block, asym, base );
                assert(sum(isnan(pbEv))==0);
            end
        end
            
        %% Tally
        function Y = tally(X)
            %tally Compute ratio after tallying across trials
            %   Detailed explanation goes here
            tmp = [{X(:,1); X(:,2); X(:,3)}];
            Y = calculate(tmp{:});
            
            function [ out ] = calculate( mEv, N, block )
                out = nan(length(mEv),1);

                for b=unique(block)'
                    mask = block==b;

                    E = mEv(mask);
                    NN = N(mask);

                    pwght = cumsum(E.*NN)./cumsum(NN);
                    out(mask) = [0.5; pwght(1:end-1)];
                end
            end
        end
        
        %% Average
        function Y = avg(X)
            tmp = [{X(:,1); X(:,2); X(:,3)}];
            Y = calculate(tmp{:});
            
            function [ out ] = calculate( mEv, ~, block )
                %avg Summary of this function goes here
                %   Detailed explanation goes here
                out = nan(length(mEv),1);

                for b=unique(block)'
                    mask = block==b;

                    E = mEv(mask);

                    pwght = cumsum(E)./[1:length(E)]';
                    out(mask) = [0.5; pwght(1:end-1)];
                end
            end    
        end
        
        %% Difference
        function Y = diff(w, X)
            tmp = [{X(:,1); X(:,2); X(:,3)}; num2cell(w')];
            Y = calculate(tmp{:});
            
            function [ out ] = calculate( mEv, N, block, w )
                %diff Summary of this function goes here
                %   Detailed explanation goes here
                out = nan(length(mEv),1);

                nH = mEv.*N;
                nT = N - nH;
                D = nH - nT;

                for b=unique(block)'
                    mask = block==b;

                    d = D(mask);

                    pwght = cumsum(d)./[1:length(d)]';
                    out(mask) = [0.5; pwght(1:end-1)];
                end

                % Squash with sigmoid
                out = 1./(1+exp(-w*out));
            end
        end
        
    end
end

