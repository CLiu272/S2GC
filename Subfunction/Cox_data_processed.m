function [ Processed ] = Cox_data_processed( X,y,varargin )
%COX_DATA_PROCESSED

 %   [...] = cox_preprocess(X,Y,'PARAM1',VALUE1,'PARAM2',VALUE2,...) specifies
    %   additional parameter name/value pairs chosen from the following:
    %
    %      Name          Value
    %      'baseline'    The X values at which the baseline hazard is to be
    %                    computed.  Default is mean(X), so the hazard at X is
    %                    h(t)*exp((X-mean(X))*B).  Enter 0 to compute the
    %                    baseline relative to 0, so the hazard at X is
    %                    h(t)*exp(X*B).
    %      'censoring'   A boolean array of the same size as Y that is 1 for
    %                    observations that are right-censored and 0 for
    %                    observations that are observed exactly.  Default is
    %                    all observations observed exactly.
    %      'frequency'   An array of the same size as Y containing non-negative
    %                    integer counts.  The jth element of this vector
    %                    gives the number of times the jth element of Y and
    %                    the jth row of X were observed.  Default is 1
    %                    observation per row of X and Y.
    %      'init'        A vector containing initial values for the estimated
    %                    coefficients B.
    %      'options'     A structure specifying control parameters for the
    %                    iterative algorithm used to estimate B.  This argument
    %                    can be created by a call to STATSET.  For parameter
    %                    names and default values, type STATSET('coxphfit').


    narginchk(2,inf);
    % Check the required data arguments
    if ndims(X)>2 || ~isreal(X)
        error(message('stats:coxphfit:BadX'));
    end
    if ~isvector(y) || ~isreal(y)
        error(message('stats:coxphfit:BadY'));
    end

    % Process the optional arguments
    okargs =   {'baseline' 'censoring' 'frequency' 'init' 'options'};
    defaults = {[]         []          []          []     []};
    [baseX cens freq init options] = internal.stats.parseArgs(okargs,defaults,varargin{:});

    if ~isempty(cens) && (~isvector(cens) || ~all(ismember(cens,0:1)))
        error(message('stats:coxphfit:BadCensoring'));
    end
    if ~isempty(freq) && (~isvector(freq) || ~isreal(freq) || any(freq<0))
        error(message('stats:coxphfit:BadFrequency'));
    end
    if ~isempty(baseX) && ~(isnumeric(baseX) && (isscalar(baseX) || ...
                                 (isvector(baseX) && length(baseX)==size(X,2))))
        error(message('stats:coxphfit:BadBaseline'));
    elseif isscalar(baseX)
        baseX = repmat(baseX,1,size(X,2));
    end


    % Sort by increasing time
    [sorty,idx] = sort(y);
    X = X(idx,:);
    [n,p] = size(X);
    if isempty(cens)
        cens = false(n,1);
    else
        cens = cens(idx);
    end
    if isempty(freq)
        freq = ones(n,1);
    else
        freq = freq(idx);
    end

    % Determine the observations at risk at each time
    [~,atrisk] = ismember(sorty,flipud(sorty), 'legacy');
    atrisk = length(sorty) + 1 - atrisk;     % "atrisk" used in nested function
    tied = diff(sorty) == 0;
    tied = [false;tied] | [tied;false];      % "tied" used in nested function

    Processed = struct('X',X,'freq',freq,'cens',cens,'atrisk',atrisk,'y',y,'sorty',sorty);

end

