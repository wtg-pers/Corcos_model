function defaultopt = nsgaopt(geom,test)
% Function: defaultopt = nsgaopt()
% Description: Create NSGA-II default options structure.
% Syntax:  opt = nsgaopt()
%         LSSSSWC, NWPU
%   Revision: 1.3  Data: 2011-07-13
%*************************************************************************


lb(1:geom.NCP)              = 0.005*ones(1,geom.NCP);
lb(geom.NCP+1:2*geom.NCP)   = zeros(1,geom.NCP);
ub(1:geom.NCP)              = 0.05*ones(1,geom.NCP);
ub(geom.NCP+1:2*geom.NCP)   = 45*ones(1,geom.NCP);

defaultopt = struct(...
... % Optimization model
    'popsize', test.population,...           % population size
    'maxGen', test.generation,...           % maximum generation
    'numVar', 2*geom.NCP,...             % number of design variables
    'numObj', 2,...             % number of objectives
    'numCons', 0,...            % number of constraints
    'lb', lb,...                % lower bound of design variables [1:numVar]
    'ub', ub,...                % upper bound of design variables [1:numVar]
    'vartype', ones(1,2*geom.NCP),...           % variable data type [1:numVar]��1=real, 2=integer
    'objfun', @objfun,...       % objective function
... % Optimization model components' name
    'nameObj',{{}},...
    'nameVar',{{}},...
    'nameCons',{{}},...
... % Initialization and output
    'initfun', {{@initpop}},...         % population initialization function (use random number as default)
    'outputfuns',{{@output2file}},...   % output function
    'outputfile', 'populations.txt',... % output file name
    'outputInterval', 1,...             % interval of output
    'plotInterval', 5,...               % interval between two call of "plotnsga".
... % Genetic algorithm operators
    'crossover', {{'intermediate', 1.2}},...         % crossover operator (Ratio=1.2)
    'mutation', {{'gaussian',0.1, 0.5}},...          % mutation operator (scale=0.1, shrink=0.5)
    'crossoverFraction', 'auto', ...                 % crossover fraction of variables of an individual
    'mutationFraction', 'auto',...                   % mutation fraction of variables of an individual
... % Algorithm parameters
    'useParallel', 'yes',...                         % compute objective function of a population in parallel. {'yes','no'}
    'poolsize', 4,...                                % number of workers use by parallel computation, 0 = auto select.
... % R-NSGA-II parameters
    'refPoints', [],...                              % Reference point(s) used to specify preference. Each row is a reference point.
    'refWeight', [],...                              % weight factor used in the calculation of Euclidean distance
    'refUseNormDistance', 'front',...                % use normalized Euclidean distance by maximum and minumum objectives possiable. {'front','ever','no'}
    'refEpsilon', 0.001 ...                          % parameter used in epsilon-based selection strategy
);



