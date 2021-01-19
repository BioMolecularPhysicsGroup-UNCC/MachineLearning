function version = splocToolsetVersion()
% write to splocLog and to screen a message at start of SPLOC
% SPLOC = Supervised Projective Learning with Orthogonal Completeness
%
% This function comments all functions that are used in the SPLOC toolset
% It also provides development notes by the original, Dr. Jacobs. 
% SPLOC toolset available at:
% https://github.com/BioMolecularPhysicsGroup-UNCC
%%                                                       development notes
% version2018a
% The initial algorithm for SPLOC was designed and developed by 
% Prof. Donald Jacobs at UNC Charlotte in the beginning of Sept 2018.
% Over the 2018 fall semester, the algorithm was refined several times, 
% and the initial SPLOC toolset was programmed by Prof. Donald Jacobs. All
% testing was done using synthetic data involving 58 dimensions describing
% 29 atoms in a 2 dimensional molecule. Alpha testing of SPLOC was 
% available in January 2019.
%
% version2019a
% Around March 15 Dr. Jacobs implemented a different version of how to 
% optimize the cost function within sploc.m. On April 14, a new method was
% completed and was suppose to be much faster and not dependent on the 
% covariance matrix generator method. At the end, the new method seems to
% have slight advantages, but nothing hands down conclusive. It may have
% some disadvanges too. The biggest advantage for the new method is that
% it has a lot of room for improvement. It seems it is ~ as good as the 
% original (hopefully a tad bit better??) but it is much more flexible and
% has not been optimized. The same name sploc.m is used in the toolset as
% a replacement, but the version2018a version of sploc.m should still have
% good lifetime of use remaining. Remaining tasks would be to assign a 
% version number to sploc.m and to fix its verbosity control (did repair 
% this due to lack of time) and, to use the new version and optimize the 
% cost function. Many possibilities are open!  
%
% version2019b
% Around October 10, Dr. Jacobs modified how the probability calculations
% for classification is done. Basic formulas for binomial distribution is
% used to provide a more absolute and correct set of probabilities. Also
% with m-modes, these probabilities multiply, and a geometric mean over
% all modes is applied. These results are context dependent and relatively
% small probabilities have meaning in an absolute sense. All necessary 
% updates were made for tables and plotting ranks. 
%
% A major update was made in the sploc function. The concept of a learning
% rectifer was added so now the projection pursuit program is optimizing 
% an efficacy score based on learning two things. 1) How to maximize the
% objective scoring function, and, 2) do this with high quality. In the
% process of learing using a leaky linear learning rectifier, the slopes
% of the learning lines are controled by quality factors. This adaptive 
% learning process through the updating of quality factors is what guides
% the projection pursuit process. 
%
% The quality factors are designed to be symmetrical between swapping of
% labels (F vs. N). In addition, the quality factors minimizes overlap 
% between a functional and nonfunctional system at the same time among all
% F & N pairs. This is handeled by monitoring statistical consistency and
% working directly with features that include mean and standard deviation.
% Cluster properties determined by geometrical consisideration along the
% projected direction are used to calculate the quality factors. Also, the
% way in which different pairs of vectors are selected to be optimized 
% within the generalized Jacobi method has been updated to include a 
% probablistic mechanism that reduces the probability for selection when 
% vector pairs converge in their relative orientation. The most active 
% vector pairs are selected with higher probability. Another way to say it
% is each distinct vector pair is placed in a queue and the longer they
% are waiting, their probability to get selected again increases as the
% process continues. Hence, the least active vector pairs are put at the 
% end of the queue but they will eventually be rechecked. This helps get
% much more consistency and higher accuracy then before, which was using 
% a greedy mechanism. In additon, the higher ranking scores have greater
% probability to increase their scores further, conistent with the the
% idea that the rich get richer. The convergence criteria is now tied to
% the consensus voting. It is important to have the proposed basis vector
% set reach equilibrium when their is statistical significance. To prevent
% cases where there is excessive crawling toward a converged solution, a 
% bail out mechanism that monitors rate of improvement in efficacy score
% was also added. Finally, the growth rate mechanism that recovers from 
% the damping mechanism is also adaptively adjusted to prevent excessive
% skipping which prevents the program from reaching a point of crawing.
%
% The input to sploc was modified substantially. First, it is now possible
% to specify the type of projection pursuit being used. For options are
% possible: 1 => maximize discriminant modes & minimize indifference modes
%           0 => maximize discriminant & indifference modes simultaneously
%          -1 => maximize indifference modes & minimize discriminant modes
%           2 => First do 1 and then do 0. This is to seed the 0-run.
%           3 => Sequencially do -1, then 1 & then do 0 for double seeding
%       
% Note: if (o == 0) && (SBV == 0) the program will not run.
% 
% It is necessary to specify an initial complete basis set as a starting 
% guess. If 0 is specified, an initial guess will be done automatically.
% Otherwise, the user can specify an initial basis set to start from. 
% This feature is important so that one can build better predictions (by 
% adding training) based on what was previously done, and one need not 
% start from scratch each time. As an additional goal, to obtain results 
% that are insensitive to initial conditions for the case where both (d&i)
% subspaces are being optimized, sploc() internally runs the projection 
% pursuit with option +1 first (to maximize discriminant modes) to get the
% optimal seeding as the initial condition for a subsequent run that is 
% performed using the mode to maximize (d&i)-subspaces. This means that an
% initial guess will get modified. It was debated whether an initial guess
% should be preserved. It is possible to bypass the option 2 by using the
% pursuit type option 0, and provide an initial guess. 
%
% The default initial guess is to consider PCA-eigenvectors based on the
% covariance matrix of the functional system, non-functional system and 
% the pooled system using a balanced average over the nonfunctional and 
% functional statistical matrices. Among the three cases, the case that 
% gives the best results (before starting the generalized Jacobi method)
% is the one that is selected to serve as an initial basis set. The new
% scoring functions are symmetrical w.r.t. to labeling, so this three-way
% check for initial conditions is a large time saver. Overall, the results
% are now insensitive to which label is used for F or N.
%
% Now mcsploc() and sploc() are both operational, fully consistent
% with respect to all other functions in the toolbox. However, mcsploc()
% maintains the original ideas of sploc(). 
% 
% version2020a
% Between December 20, 2019 to Feb 1, 2020, clean up work was done to
% organize the code better for human understanding, adding more comments
% that explain rationale for processes/methods, and to optimize it a tad.
% and to prepare the code for public distribution at the BMPG github site.
% Before this time, 99% of the effort was to write code, not optimization.
% In this period of time, algorithms were parameterized more optimally by
% incorporating adaptive feedback loops wherever possible, existing and/or
% new algorthms were improved/created whenever opportunity was found. This
% work created more work in varifying that functionality is maintained 
% across all the functions in the toolset. Therefore, additional test runs
% were made to look for bugs and to fix them, which also highlighted the
% areas where optimization needed further improvement. In addition, a test
% code for multi-class version of SPLOC was fully implemented, but the
% mcsploc() function is still under development. The main purpose of doing
% this was to test the generalities of the methods/algorithms being used 
% in sploc. All methods in sploc now generalize to the multi-class case. 
% The upshot of all the technical work done in this period of time is that
% the code is NOT considerably faster, but SHOULD be more accurate with a
% a small speed increase while the results SHOULD be more consistent. As a
% BONUS it has been shown that multiclass discrimination via SPLOC is not
% only feasible, but efficient. However, a battery of tests is required by
% BMPG members before that part of the code can be released to the public.
%
% To make sploc functions generalize to multi-class discrimination, the 
% output format had to be generalized, so this required a complete change
% in datastructures. For binary classification, the data structures that
% are employed now are arguably more complicated than need be. However, 
% this was intentionally done so that all functions in the SPLOC toolset
% immediately generalize to the problem of multi-class discrimination. The
% coloring scheme and details of how the results are presented have been
% upgraded based on feedback from BMPG users with binary classification,
% and preliminary results from the multiclassification results using the 
% mcsploc() function and test datasets.  
%
% The main reason/impetus that started this bout of update was to improve
% the quality factor definition. It was noticed in the 2019b version that
% a domain user could select discriminate modes that exceed both selection
% and consusus power thesholds, but nevertheless lead to poor clustering
% (no separation in feature space) and actually reduced the quality of the
% classification probability predictions (soft classification) on unknown
% systems. To work around this lenient criteria, the user had to visually
% look at the results of all discriminate modes (or d-modes) for many runs
% and manually remove poor performing d-modes. This begged the question
% that certain linear combinations of the d-modes (good & bad ones) could
% still be better than the remaining good modes a user selects. All this
% fuss then diminished consistency, since it was then necessary to run
% dozens of cases, manually checking everything to determine the best
% cases. To confirm, a user had to repeat this same process at least a few
% times to make a strong conclusion. With scripts this could take a couple
% of days of work for a domain expert when working with molecular dynamics
% data that would take a half-a year to collect. Still, not bad, but, this
% is arguably a messy process simply because the quality factor is missing
% information about the quality of the clustering itself. The human then
% has to make inference based on looking at the results and comparing. As
% such, the criteria for a projection to be a d-mode has been extended to
% make clustering in feature space 100% separable as a guarantee. As an 
% upshot, the need to do extensive post comparative analysis has been 
% greatly reduced because sploc builds this requirement in its exploration
% as part of projection pursuit. Thus, one sploc result is doing 99.5% of
% all the work. It is still recommended a domain user perform extensive 
% post analysis, since the bar of what is good/bad has shifted upward.
% How did this time saving occur? 
%
% Now the quality factor is ZERO unless there is 100% separation between
% all functional and nonfunctional training points within feature space.
% For example, if there are 10^9 F-points and 10^9 N-points, and 1 F-point
% is infinitesimally not separable from one N-point, the quality is 0.
% Even if you were comparing 10^1000 functional points, and 10^1000
% nonfunctional points, and all were separated except 1 functional did not
% separate, quality will be ZERO. Yes, this is very restrictive, but it is
% even MUCH MORE servere than this! Since feature space deals with means
% and standard deviations of packets of information, uncertainties are 
% placed on the mean and standard deviation, depending on the sample size
% of the packet. A conservative estimate of a BAD fluctuation is added
% so that there actually needs to be some "space" between clusters that 
% differentiate two different things. That space can be smaller as the 
% uncertainty in the moments decrease. This severity tremendously IMPROVES
% results when found, but the problem then happens that most of the time
% a solution is hard to find a projection that yeilds full separation. To
% overcome this problem, the concept of nucleation from physics was added
% as part of the optimization process. Clusters are hard to form at first
% due to a nucleation barrier. From chemistry, it is well known that to
% get over a nucleation barrier one needs a catalysis. This is done by
% allowing for precursory findings to influence the search. If clustering
% seems to improve (but not fully consistent) the projections will be 
% optimized, to help the learning process to overcome the threshold, but, 
% the projection will not be declared as d-mode until the hard threshold 
% is satisified. Said another way, it is like how a professional sports 
% team works. Pick only the best of the best, but help train the not best
% until they become one of the best of the best. This is what it means to
% learn, and is essentially the essense of the meaning of SPLOC. So not 
% finding awesome solutions is no longer an issue due to the prenucleation
% tricks added to the code. The upshot of this modification is that when
% no discriminate modes are found, it indicates that there is at least one
% system that is NOT different from a set of systems that are claimed to
% be different. Obviously the statement is wrong, because to prove it is 
% wrong only requires ONE counter example. Hence, with only ONE exception,
% among potentially a googol worth of tests, makes that statement wrong. 
% For a user to make this statement with reasonably certainty, requires
% the user to increase statistics. The REASON for making the condition so
% difficult (based on statistics alone) is to enforce the general saying
% of Dr. Jacobs that says: "STATISTICS DO NOT LIE: PEOPLE DO". If a user
% is to say that two systems are different based on a feature, it should
% be a true statement for at least the TRAINING set!!! 
%
% A summary of the criteria that makes a projection a discriminate mode:
%
%  If: quality   > minimum quality threshold   => clear separable clusters
% and: consensus > minimum consensus threshold => statistical significance
% and: selection > minimum selection threshold  => exceeds signal to noise
%
% then a projection is declared as a discriminate mode based on the 
% sampled data on hand, which determines the uncertainties. Even if all
% thresholds are met, does not mean one has a true model, but it is at 
% least consistent with the data on hand.
%
% Plotting: Now six graphs pop out of sploc.
% GRAPH 1:
% selection power,
% consensus power,
% quality factor,   
% -------------------> these three plots on 1 graph help a domain user see
% the usefullness of a projection. Discriminant and indifference modes are 
% now colored in red and blue respectively, and yellow is still used to
% indicate undetermined. A domain user can see which of the properties are
% or are not satisified.
% GRAPH 2:
% efficacy,
% within similarity for indifference modes
% within similarity for discriminant modes
% ------------------------------------------> provides nice summary for
% how good projections are. Efficacy could be used as a sorting index, but
% selection power is still selected as the sort ordering indexing scheme 
% for sploc projections. The exact sorting method is more complex for the 
% multi-class case, but is designed to be fully consistent with the binary
% classification case, for which there is no ambiguity. 
% 
% other optimization tricks.
% Better use of time delays and importance sampling has been added. These
% can be tracked when sploc is run with verbosity of 3 to help debug and
% understand what the code is doing. Essentially, spining pairs of vectors
% in their plane is needed to find better solutions, but the trick is to
% consider pairs that will yield better results more consistently than 
% other pairs that do not yeild improvements. By tracking information on
% last time a pair of vectors were spun, and the gains (if any) that were
% made on a previous spin, and which vector tends to gain more than other
% vectors, it is possible to increase overall sucess rate. However, to do
% this, one must skip spinning some vector pairs as well, and waiting to 
% spin them at a later time. All this is done by forming a queue, and
% selecting certain reference vectors with more frequency than others, and
% the conjugate pairs of those vectors with more frequency than others, 
% all based on statistical results of previous spins. This developes a 
% much faster algorithm, but does have some chance of biasing the result.
% In fact, the biasing is needed to make the process run faster, since 
% checking all pairs all the time leads to a slow program. To mitigate 
% against missing pairs that were incorrectly identified as uneffective, 
% there are a few mechanisms that look for the "diamond in the rough".
% This however, increases the total time of the algorithm, but it also 
% decreases the variance in solution quality. All the paramters, which 
% are almost always adaptive based on data-driven sampling, have been 
% set based on heristic ideas that come down to how a system learns, 
% whether a human learning technical information, a skill, a sport or AI.
%
% In addition to tweaking the importance sampling monitoring processes,
% and the nucleation trick, another random mechanism was added, to help
% shake the system away from getting no solution, or to help push out of
% a local minimum in order to find a better solution compared to the one
% on hand. This is tricky, because shaking the system out of equilibrium
% can lead to a long process of the system relaxing back to equilibrium,
% which occurs when the shaking is too great so that convergence criteria
% cannot become satisified. This is essentially saying the perturbation on
% a dynamic system is too frequent and/or has too great of an amplitude 
% compared to the intrinsic relaxation time of the system. The shaking of
% a system away from its equilibrium in the form of a stochastic process,
% analogous to a Langevin approach, has an advantage of strengthening
% the ability of a solution to hold up to improvisational performance by 
% allowing the system to explore more projections than it would have
% otherwise. All this is implemented by randomly rotating projections in
% the undetermined subspace only when u-modes have a low precursor 
% indication of nucleation. Note that by definition the efficacy of a 
% u-mode is zero, and as such, the u-subspace in terms of efficacy is
% fully degenerate. The number of such vectors that are randomized, and 
% the amplitude of rotation (mean deviation in angle between initial and
% final) are all adaptively adjusted based on a data driven scheme. A 
% conservative approach was taken here because it was found that often, 
% not using the random mechanism performed better. However, this is not 
% always true. If you run sploc say 10 times, you might get a better and
% faster solution 3 times out of 10 if the random rotation mechanism is 
% turned off. But with the adaptive format now implemented, it is likely 
% 10 good solutions will be found out of 10 tries, but 3 of them might be
% as good as the 3 that might have been found without it. On the other 
% hand, the 7 other solutions without random rotations could be relatively
% poor solutions compared to the best. In general, the random rotation 
% mechanism helps much much more than it deters. One way to think about 
% what the random rotations in the u-space is doing is it is like giving
% bateria some stress where it is likely that most will die, but survivors
% will be stronger. The reason why the random rotations make the solutions
% stronger, is because has d-modes or i-modes spin with randomized u-modes
% they are spinning in planes that were not spun before, and hence there
% is more exploring for vector space. Thus randomization is good, if the
% frequency is not too high, and for moderate to low amplitudes. In order
% to reach a good solution, despite the presense of the stochastic noise
% term that pushes the system away from equilbrium, the noise term is 
% turned off on the end-game of the optimization process. The noise is
% strongest in the beginning of the process and decreases as one settles 
% to a good solution (the meaning of good solution is in a global sense).
% As such, it is also valid to interpret the random rotations as a bona 
% fide diffusive process (or Brownian motion) that has an amplitude that 
% decreases as a solution coverges in the same way that is applied using
% simulated annealing. The random fluctuations that are imposed on the 
% undetermined subspace is large when nothing is known about the system,
% and shrinks as more is known. Furthermore, if the success rate is high
% using the generalized Jacobi method, the random rotation is not really 
% needed, so it is slowed down. In constrast, when the success rate of a
% vector pair spin is low based on spinning pairs of vectors, there is not
% much to loose because the system is literally spinning its "wheels" not
% going anywhere. As such, the frequency and amplitude of random rotations
% are tied to convergence criteria and success rates, and quality of the 
% u-modes that are being spun.
%
% Note that the random rotations are implemented using the Cayley method.
% The number of rotations made, based on amplitude and frequency should be
% interpreted as a diffusive process, and the results of such a process is
% reported in that language. 
% 
% The getBasisVecSpectrum.m function was updated so it can handle any 
% number of modes, and it will also provide an indexing map between an
% input set of basis vectors, to the sploc modes. These features were 
% desirable to make comparisons of other machine learning approaches 
% against sploc, and, in fact, this allows those other methods to have
% an context dependent evaluation on them. Also, as mentioned above, a
% completely new data structure was developed to handle eventual multi-
% classification applications using mcsploc().
% 
% version2020b
% Modifications were manly made in sploc() initially motivated by a bug
% that Chris Avery spotted. A patch was first made to remove a potential
% infinite loop, which was thought to be responsible for the occasional
% unexplained excessive slowness that sploc() has. Searching through the
% code further revealed a few more bugs that were fixed while fine tunning
% was applied in the optimization algorithms. Much of the code was then
% streamlined, and a [gremlin] was expelled. This was another place that
% sploc() could get excessively slow. Most important changes are the 
% introduction of mode-type probabilities, initial probabilities are
% randomized, sequential resets are made in a simpler way, including 
% the stopping of anemic performance, and the gremlin was found related
% to skipping rotations. By adding a new concept of quorum the gremlin was
% expelled. a LOT of MINOR & MAJOR changes were made [~1200 lines] keeping
% the basic ideas the same, but adding several NEW features that help the
% importance sampling. The front end of sploc() was modified slighly too.
%
% Current options:
%
%      splocResults = sploc(o,X,baseFname,trait1,trait0,verbosity,vT)
%                           ^ ^---> no initial guess for basis set
%                           |-----> SPLOC pursuit objective
%
% X = 0,1,U0  0 => internally the best PCA pooling case will be used.
%             1 => identity matrix will be used for intial guess. 
%            U0 => user defined basis set can be used as an initial guess.
%
%                                  AVAILABLE MODES
%       o = 2 => maximize discriminant modes & MINIMIZE indifference modes
%                RANU is modified to give inversion of i-modes
%                sets reference value for sVS = 0.1 = set_sVS (standard)
%       o = 1 => focus on maximizing discriminant modes
%                sets the reference value for sVS = 0.01 = set_sVS
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Dd
%NORMAL o = 0 => maximize discriminant & indifference modes simultaneously
%                sets reference value for sVS = 0.1 = set_sVS (standard)
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Dd
%       o= -1 => focus on maximizing indifference modes
%                sets the reference value for sVS = 1.0 = set_sVS
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Di
%       o= -2 => maximize indifference modes & MINIMIZE discriminant modes
%                RANU is modified to give inversion of d-modes
%                sets reference value for sVS = 0.1 = set_sVS (standard)
% 
% RARE^ use o = 1 or -1 when a bias toward d-modes or i-modes is desired.
%               When the bias for d-modes increases, the bias for i-modes
%               is unaltered. Thus, i-modes only get removed because a 
%               d-mode can be created. Vice versa, when the bias for an
%               i-mode increases, the bias for d-modes is unaltered. Thus
%               i-modes only get removed because a i-mode can be created.
%
% RARE^ use o = 2 when a bias for d-modes and AGAINST i-modes is desired.
%               Generating i-modes decresses efficacy.
%
% RARE^ use o = -2 when a bias for i-modes and AGAINST d-modes is desired.
%               Generating i-modes decresses efficacy.
% 
% Outside of these small changes on the input generalizations to sploc()
% the rest the changes made in this last release is backward compatible.
%
% The most important added concept was instead of trying to get all the 
% modes to converge accurately at the same time, the program focuses down
% on the side of spectrum MOST important to converge first, and the rest
% of the modes follow like a domino effect. A relative weighting between 
% d-modes and i-modes was incorporated which was critical in finding
% d-modes and i-modes when desired. This concept is implemented as an 
% adaptive multiplicaive factor on the relative weighting between the
% d-modes and i-modes. This was always a source of confusion, but now the
% program internally adaptively optimizes this factor. In additon, the 
% selection process of the first mode of a pair, called j1, was further
% optimized by targeting parts of the mode spectrum of interest. This is
% implemented by selecting modes in order of preference according to their
% rank order dictated by the scoreing function (interesting index). Also  
% an adaptive weighting is used to blend random ordering with rank sort
% ordering. Additonal optimization was done in terms of how successful
% rotations are defined. Now, efficacy must be increased above a minimum
% level for success to be defined. Additional modifications were made in 
% how resets are counted after the restructuring was made, and how big 
% nTimes is set. In additon, the classifer probabilities were changed 
% to a root mean square average plus a small modification of these values
% using a weigthted geometric mean using a hybrid model. At the end, sploc 
% has much faster convergence, runs faster and achieves greater accuracy. 
% Moreover, sensitively to initial conditions has been greatly dimenished.
% The 2020b release was finished on May 2, 2020.
% 
% version2020c
% "... and the beat goes on".
%%                                                     directory structure
% Toolset:
%
% I/O directory structure:
% 
% current directory --> sub-directories: splocLibrary
%                                        input
%                                        splocLog
%                                        training
%                                        basisComparison
%                                        classification
%                                        analysis
%                 
%%                                                       Test code scripts 
%
% script testDriver4control
%    --> Test functions related to controling the sploc functions 
%
% script testDriver4input
%    --> Test functions related to the input sub-directory. 
%
% script testDriver4training
%    --> Test functions related to the training sub-directory.
%
% script testDriver4basisComparison
%    --> Test functions related to the basisComparison sub-directory.
%
% script testDriver4classification
%    --> Test functions related to the classification sub-directory.
%
% script testDriver4analysis
%    --> Test functions related to the analysis sub-directory.
%
% script testDriver4workflow01
%    --> Test SPLOC in a scenario that emulates a protein design workflow
%
% script testDriver4workflow02
%    --> Test MCSPLOC in a scenario that emulates protein design workflow
%
% script testDriver4workflow01i
%    --> same as testDriver4workflow01 except a swap is made in functions 
%        getMultivariateStats4sploc()  --->  getBoostedMVStats4sploc() 
%%                                             I/O characteristics summary
% INPUT: 
% Many types of files can be read by a variety of sploc functions. Each
% such file must be formated according to specifications. In many cases, 
% the output files generated by some sploc functions will serve as input
% files for other sploc functions, which will usually have dependence. For
% example, under supervised training on various datastreams that are known
% to be functional or nonfunctional, a basis set of vectors is optimized
% to successfully classify an unknown datastream. At a later time, new
% acquired data can be immediately classified by reading in an input file
% that contains the relevant basis vectors for discrimination, without 
% performing new training. Hence, the output from training becomes the
% input for classification along with more data. 
%
% OUTPUT
% verbosity: Specifies amount of intermediate processing steps to report.
%            Functions exist to output certain data, while verbosity has
%            no connection to that type of output and how it is reported.
%       ---> verbosity extracts hidden information within SPLOC functions.
%       ---> verbosity controls optional intermediate results ONLY:
% default: 0 => minimal checklist report in sploc log file.
% process: 1 => same as 0, plus detailed report in specialized log file.
% summary: 2 => same as 1, plus writing figures to show key relationships.
% display: 3 => same as 2, with figures displayed on the screen, sometimes
%               with a pause, and/or with additional output printed to the 
%               command window. Useful to understand how the code works!
%       ---> All printed figures have the same file type (fig, png, etc).
%       ---> Specialized log files are written in separate directories.
% ------------------------------------------------------------------------
%%                                                       control functions
%
% function initializeSPLOC(gFileType,verbosity)
%      --> sets all global parameters used by functions shared across the
%      --> SPLOC toolset. The splocLibrary is accessed. This function
%      --> should always be the first sploc function called. 
%
% function longFileName = getOutputFileName(subFolder,fName)
%      --> creates subfolder and/or augments folder as a prefix to the 
%      --> output file name. 
%
% function inFileName = getInputFileName(targetFolder,fName)
%      --> makes targetfolder a prefix to the input filename (fName) as it 
%          checks if the path/filename exists, leading to an error if not. 
% 
% function version = splocToolsetVersion()
%      --> generates the version of the MATLAB SPLOC toolbox 
%
% function verbosity = setVerbosity(verbosity)
%      --> enforces verbosity input to have admissible values {0,1,2,3} 
%
% function setGraphicFileType(gFileType)
%      --> resets the graphics file type after initialization
%
% function matrixFormat = setDataMatrixFormat(pType,dim)
%      --> defines the packing characteristics for the vector space 
%
%%                                                       utility functions
% Remarks: Due to their ubiquity, a test driver is not needed for utility
% functions. Using a utility function does not get recorded in log files. 
%
% function fDate = getDateString()
%      --> makes a string out of the current date. Used for splocLog file.
%
% function lineName = dividerLine(str)
%      --> creates a string variable for dividing information within a 
%      --> file. Extensively used in log files.
%
% function rmsip = getRMSIP(U,V)
%      --> calculation of RMSIP between basis sets U and V for subspaces.
%
% function rmsipSequence = getRMSIPsequence(U,V)
%      --> yields a sequence of root mean squared inner products as DIM(U)
%      --> is scanned, while the subspaced defined by V is held fixed. 
%
% function reducedC = reduceCmatrix(C,pType,dim)
%      --> barebones utility to reduce the C-matrix that has vectorized
%      --> components to a smaller matrix that uses inner products to 
%      --> create a scalar description of the same data, and not vectored.
%
% function U = getDiscriminantSBV(splocResults)
%      --> The selection basis vectors spanning the discriminant congruent 
%          subspace is extracted from a splocResults data structure. U is
%          either a matrix of dimension (nVariables x nModes) or a cell
%          array of such matrices that will generally have varying nModes.
%
% function U = getUndeterminedSBV(splocResults)
%      --> The selection basis vectors spanning the undetermined congruent 
%          subspace is extracted from a splocResults data structure. U is
%          either a matrix of dimension (nVariables x nModes) or a cell
%          array of such matrices that will generally have varying nModes.
%
% function U = getIndifferenceSBV(splocResults)
%      --> The selection basis vectors spanning the indifference congruent 
%          subspace is extracted from a splocResults data structure. U is
%          either a matrix of dimension (nVariables x nModes) or a cell
%          array of such matrices that will generally have varying nModes.
%
% function [modePROJ,modePDF] = generateModePDF(A,basisVectors)
%      --> Project data matrix onto basis vectors to construct PDF per 
%      --> mode. This function uses ksdensity() and can be generalized to
%      --> used PDFestimator. modePROJ is the "x" and modePDF is PDF(x).
%      --> FIX ME: An option should be included to use either the built in
%      --> ksdensity() function or the much better PDFestimator. 
%
% function colorMatrixTool(M,varargin)
%      --> a generic tool not associated with sploc used for coloring 2D
%          matrices such as covariance matrices and various heat maps. 
%          This function is integrated into many sploc functions, however.
%
% function map = lbmap(n,scheme)
%      --> Returns different color spectrums inspired by Light-Bertlein
%      --> colormap. Used as a required component within colorMatrixTool.m
%
% function overlap = getOverlap(x1,pdf1,x2,pdf2)
%      --> calculates overlap between two probability density functions
%      --> where overlap = [2 - integral |pdf2 - pdf1| *dx ]/2
%
%%                                              content of data structures
% dataMatrixInfo. <-- data structure
% dataRefName    = reference name for A matrix data with similar traits
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%          pType = packing type (format) applied to the vector space
%            dim = # of components in local vector (e.g. x,y,z for 1 atom)
%              n = number of data matrices
%           A{:} = cell array for the data matrices in the rowVar format
%    nSamples(:) = array for the number of samples in each data matrix
%
% traitData.  <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%  nVariables = number of variables defining size of vector space
%       pType = packing type (format) that is applied to the vector space
%         dim = # of components in a local vector (e.g. x,y,z for 1 atom)
%           n = nis = number of independent subsamples in this collection
%          mu = cell array for mean vectors representing the system
%          cM = cell array for cMatrix (covariance, correlation, etc) 
%  sampleSize = sample size for partitioning data matrices (a fixed value)
%     nDtotal = total # of data samples = nis*sampleSize
%
%
% splocResults. <-- data structure
% splocResults.sType       = type of spectrum: MCSPLOC
% splocResults.pursuitType = 1,0,-1  => d, d&i, i  set DEFAULT VALUE = 0
% splocResults.pType       = packing format describing the vector space
% splocResults.dim         = # of components in a local vector
% splocResults.baseFname   = base file name
% splocResults.SBV         = selection basis vectors
% splocResults.vT          = voting threshold to establish consensus
% splocResults.efficacy    = quantifies the ability for SBV to cluster
% splocResults.Dd          = # of discriminant only modes
% splocResults.Ddi         = # of (discriminant & indifference) modes
% splocResults.Du          = # of undetermined modes
% splocResults.Di          = # of indifference only modes
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% splocResults.EEVd = efficacy eigenvalues for discriminant congruence
% splocResults.EEVi = efficacy eigenvalues for indifference congruence 
% splocResults.USVd = conditional upper similar discriminant-eigenvalues
% splocResults.USVi = conditional upper similar indifference-eigenvalues
% splocResults.LSVd = conditional lower similar discriminant-eigenvalues
% splocResults.LSVi = conditional lower similar indifference-eigenvalues
% splocResults.QEVd  = conditonal quality discriminant-eigenvalues
% splocResults.QEVi  = conditonal quality indifference-eigenvalues
% splocResults.SEVd  = conditional selection discriminant-eigenvalues
% splocResults.SEVi  = conditional selection indifference-eigenvalues 
% splocResults.CEVd  = conditional consensus discriminant-eigenvalues 
% splocResults.CEVi  = conditional consensus indifference-eigenvalues 
% splocResults.Cind  = congruency indicator (2,1,0,-1)
%                    2 => discriminant and indifference congruences
%                    1 => projections for discriminant congruences
%                    0 => projections that are undetermined
%                   -1 => projections for indifference congruences
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% systemRanking. <-- data structure
%     rankType = ('dataStream','moments')
%    baseFname = base file name useful for ID and spawning more file names
%    nXsystems = number of systems assigned a classification rank
%                X represents any system (0,1,unknown) that is ranked
%                1-system => from training set labeled as functional
%                0-system => from training set labeled as nonfunctional
% dataRefNameX = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                X represents any system (0,1,unknown) that is ranked
% dataRefNameF = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                Applies to the 1-systems that were used in a training set
%                F represents 1-systems that are regarded as functional
% dataRefNameN = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                Applies to the 0-systems that were used in a training set
%                N represents 0-systems that are regarded as nonfunctional
%       nModes = # of discriminant modes used to classify unknown systems 
% ------------------------------------------------ list of cell arrays {:}
% systemXname1 = traitData.mMatrixName         for rankType =    'moments'
%              = dataMatrixInfo.dataMatrixName for rankType = 'dataStream' 
% systemXname2 = traitData.cMatrixName         for rankType =    'moments'
%                X represents any system (0,1,unknown) that is ranked
% -------------- list of arrays that reveal classification characteristics
%     pFave(j) = estimated Likelihood the j-th X-system is functional. 
%     pFstd(j) = standard deviation in the estimated likelihood.
%        wm(k) = average weight for k-th mode taken over all training data
% ---------------------------------------------------------- double arrays
%    pXF(j)(k) = MODE likelihood the j-th X-system is functional in the
%                k-th mode without comparing to known N-systems.
%    qXN(j)(k) = MODE likelihood the j-th X-system is not nonfunctional in
%                the k-th mode without comparing to known F-systems.
%    pqX(j)(k) = MODE likelihood the j-th X-systm is functional based on 
%                being similar to {F}-systems & dissimilar to {N}-systems
% ------------------------------------------------ discriminant SBV matrix
%            U = nV x nModes matrix containing nModes used to discriminate
%                data, where each mode lives in a nV dimensional space. 
%
% RMSFinfo.       <-- data structure
% dataRefName    = inherits name of the A matrix used to create the RMSF
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%              n = number of RMSF descriptions
%        rmsf{:} = cell array for the RMSF function & length(RMSF) = #DOF.
% REMARK:          the index for DOF automatically ranges from 1 to end
% 
%
% featureX. <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%   nXsystems = number of systems being projected into feature space
%               X represents any system (0,1,unknown) that is ranked
%                 1-system => from training set labeled as functional
%                 0-system => from training set labeled as nonfunctional
%      nModes = # of discriminant modes contained in U.
%   nFeatures = 2xnModes = number of distinct features
%     Fmatrix = nFeatures x nXsystems  
%%                                                                   input
%
% function [nFiles,varargout] = readFileNameList(inputFname,varargin)
%      --> Reads an input file specifying a list of file names with 
%      --> identifiers so that these other files can be subsequently read
%      --> using other read functions, but with their data classified as
%      --> either functional, nonfunctional or unknown {F,U,N}. Optionally
%      --> the file list can be read without identification. In addition,
%      --> the output from this function can be restricted to only of a
%      --> certain class.
%
% function [nFiles,pID,fileList] = readPropertySpecs(inputFname,varargin)
%      --> Similar purpose as readFileNameList. From the fileList output,
%      --> a list of file names will be obtained from the input file, and
%      --> these files specify the type of data that is to be used for 
%      --> training via specifying properties of the files through a 
%      --> property identification (pID). 
%
% function traitData = getMultivariateStats4sploc(dataMatrixInfo, ...
%                      sampleSize,verbosity,cType)
%      --> converts data matrices into covariance or correlation matrices.
%      --> User has a way to partition the sampling within a data matrix
%      --> to account for consistency in statistics by trying to remove
%      --> overfitting to outliers. This function will produce more files
%      --> in the input directory because this statistical information can
%      --> be used as necessary input to run further analyses.
%
% function traitData = getBoostedMVStats4sploc(dataMatrixInfo, ...
%                      nPairs,verbosity,cType)
%      --> converts data matrices into covariance or correlation matrices.
%      --> User has a way to partition the sampling within a data matrix
%      --> here the data is shuffled in order and split in half. Each half
%      --> is used with half of the samples. But each half gives twice the
%      --> output, which is why the user controls the number of pairs of
%      --> ouput covariance or correlation matrices. 
%      --> to account for consistency in statistics by trying to remove
%      --> overfitting to outliers. This function will produce more files
%      --> in the input directory because this statistical information can
%      --> be used as necessary input to run further analyses.
%
% function [dataMatrixInfo,T] = readDataMatrices(dataRefName, ...
%                               dataMatrixName,mFormat,fType)
%      --> read list of file names for data matrices, then read those data
%      --> matrices. Serves as a critical component in creating traitData
%      --> for other sploc functions. Output is defined by: (1) number of
%      --> vector components in the data matrix corresponding to number of
%      --> variables (or DOF). (2) Number of samples in the data matrices.
%      --> (3) The data matrices. Options allow row-column or column-row
%      --> formats to be handled with ease. Note that regardless of the 
%      --> input format, the output format of a data matrix A has its rows
%      --> represent variables and its columns represents samples. Since
%      --> all functions in SPLOC involving data matrices use the rowVar 
%      --> format convention, any necessary conversion is done on step 1. 
%
% function [traitData,T] = readTraits4sploc(dataRefName, ...
%                          cFname,mFname,nSamples,mFormat)
%      --> guided by the specification of file name lists, and character
%          of the statistical significance (nSamples) and vector space 
%          format (mFormat) the data matrices/vectors are read from files
%          and packaged into the traitData data structure that is used to
%          sploc the data. A summary table showing origin of data is given
%
% function writeSPLOCresults(fileName,splocResults)
%      --> Write splocResults (a data structure) to the file: fileName
%
% function splocResults = readSPLOCresults(fname,fullPath)
%      --> Retrieves relevant SPLOC information from specified input file.
%      --> The fullpath toggle of (0,1) allows for relative/absolute paths
%
%%                                                                training
%
% function getDefaultVoteThreshold((nS0,nS1,nDOF)
%      --> helper function to calculate the default vote threshold based
%      --> on heristics to maintain consistency in definition throughout 
%      --> the program. Note nSX = nSamples*sqrt(nX) = # of X-classes.
%      
% function splocResults = sploc(o,U0,baseFname,trait1,trait0,verbosity,vT)
%      --> Determines selection basis set associated with 3 congruency
%      --> partitions: {discriminant, undetermined, indifference}. The
%      --> method requires trait1 and trait0 to characterize similarities
%      --> within trait groups and differences between trait groups. An
%      --> initial guess for the complete basis set can be given as input.
%      --> various objective pursuit types can be specified as well. The 
%      --> sploc() function is the foundation for the SPLOCtoolset. 
%
% function mcsplocResults = mcsploc(o,U0,fName,traitL,cID,verbosity,vT) 
%      --> similar to sploc() but works for labled data with multiclasses.
%      --> The training set is defined in terms of one set of traits in a 
%      --> particular order, and corresponding array for class labels.
%
%
% function [splocResults,varargout] = ...
%        getBasisVecSpectrum(U,fname,trait1,traitORcID,vT)
% USAGE:
%
%      [splocResults,mapU2S,mapS2U] = ...
%                       getBasisVecSpectrum(U,fName,trait1,trait0,vT); 
%      --> Given a proposed basis set U: Based on a training dataset for 
%      --> supervised classification using the selection and consensus 
%      --> spectrums as well as quality factors the basis vectors are 
%      --> indexed and placed into appropriate congruent subspaces. 
%      --> Importantly, no optimization on the U basis vectors is
%      --> performed. The purpose of this function is for evaluation only.  
%      --> Note that U need not be a square matrix. Any number of modes
%      --> can be evaluated size(U) = [nDOF,nModes], where nDOF = # of 
%      --> degrees of freedom or variables, and nModes = # of modes. The
%      --> output is variable. Along with splocResults, the mapping from 
%      --> the old indexing (determined by initial input) to the new 
%      --> index is outputed, or from the new to the old index.
%      
%      mcsplocResults = getBasisVecSpectrum(U,fname,traitL,cID)
%      --> this version is essentially the same as above, except the input
%      --> has a different format. The training set is defined in terms of
%      --> one set of traits in a particular order, and there is an array
%      --> that labels each dataset into a class. This will allow for 
%      --> multiclass classication. 
% 
% function [splocResults,mapU2S,mapS2U] = ...
%          getBVStwoClasses(U,baseFname,trait1,trait0,vT)
%      --> guts of the  getBasisVecSpectrum.m  function defined above. 
%      --> This function is called within getBasisVecSpectrum() whenever
%      --> the process of binary classification is invoked. 
%
% function [splocResults,mapU2S,mapS2U] = ...
%          getBVSmultiClass(U,baseFname,traitL,cID,vT)
%      --> guts of the  getBasisVecSpectrum.m  function defined above. 
%      --> This function is called within getBasisVecSpectrum() whenever
%      --> the process of multiclass classification is invoked. However,
%      --> this program is not finalized yet.
% 
%%                                                         basisComparison
%
% function [msipMap] = subspaceComparison(cSBV1,cCIP1,verbosity,fname)
%      --> calculates the mean squared inner product matrix and plots the
%      --> results depending on verbosity for two modes of operation. 
%      --> MODE 1: 
%      --> When cSBV1,cCIP1 are specified as shown above, the objective is
%      --> to quantify the consistency of a single calculation type. Note
%      --> that the cell array of various basis vectors hold variations of
%      --> a given calculation. It is up to the user to determine what is
%      --> being compared, but often the data comes from an ensemble of 
%      --> similar solutions.
%      --> MODE 2:
%
% function [msipMap] = subspaceComparison(cSBV1,cCIP1, ...
%                                         cSBV2,cCIP2,verbosity,fname)
%       --> Same as above, but two different sets of basis vectors can be
%       --> compared. This allows one to compare variations within two
%       --> types of calculations at the same time while comparing two
%       --> different types. 
%%                                                          classification
%
% function [ranking,T] = classifyMoments(fname,traitU,traitF,traitN, ...
%                                        U,verbosity)
%       --> Calculates PROB[unclassified system is functional] based on 
%       --> the average and standard deviation of the data projected onto
%       --> mode directions. It is assumed that the probability density 
%       --> functions for the data are described by a Gaussian probability 
%       --> density. Overlap integrals are used together with activation 
%       --> functions to calculate the probability that a system viewed 
%       --> along a particular basis vector direction has functional 
%       --> properties. A weighted average over these probabilities per 
%       --> basis vector is taken where the weights depend on separation 
%       --> between known functional and nonfunctional systems.
%
% function [ranking,T] = classifyMoments(fname,traitU,traitF,traitN, ...
%                                        U,verbosity)
%       --> Calculates PROB[unclassified system is functional] based on an 
%       --> estimated probability density function using either ksdensity
%       --> or PDFestimator function. Overlap integrals are used together
%       --> with activation functions to calculate the probability that a 
%       --> system viewed along a basis vector direction has functional 
%       --> properties. A weighted average over these probabilities per 
%       --> basis vector is taken where the weights depend on separation 
%       --> between known functional and nonfunctional systems.
%
% function featureX = getFeatureVectors(traitX,U)
%       --> Calculates feature vector for each system defined in traitX
%       --> There are two features per basis vector given in U per input 
%       --> system X. The features are the mean and standard deviation per
%       --> selection basis vector. A system of M dimensions gets mapped 
%       --> into feature space of 2M dimensions. However, the dimension of
%       --> the discriminant subspace, within feature space will generally
%       --> be much smaller than number of DOF (nDOF) in the system. 
%%                                                                analysis
%
% function plotCongruencySpectrum(splocResults,verbosity)
%      --> plot the selection power and consensus power summary bar graphs
%
% function setPlotCongruencySpectrum(k1,k2)
%      --> k1 = (1,-1) for plotting information upright, or top/bottom.
%      --> k2 = (0,1) for (not including, including) geometric mean line.
%
% function dataOutput = dataStreamProjection(fname,dataInput,U)
%      --> Projects A-matrix data into a subspace defined by basis vectors
%      --> in U. The output format is the same as the input format, which
%      --> is defined by the  dataMatrixInfo  data-structure. Note that
%      --> the output A-matrix will generally be much shorter than the 
%      --> original A-matrix because we are using a subspace defined by U.
%      --> In particular the basis vector components are essentially the
%      --> same as principal components or loadings. 
%
% function dataOutput = dataStreamFilter(fname,dataInput,U)
%      --> Filters out information from the A-matrix that is outside the U
%      --> subspace. A projection operator is created as U*U' which is 
%      --> applied to A. The math   Afiltered = (U*U')*Ainput   shows that
%      --> the outputted filtered A-matrix has the same dimensionality as
%      --> the original A-matrix. Notice the output format in dataOutput 
%      --> is the same as the input format, which is defined by the 
%      --> dataMatrixInfo  data-structure. 
% 
% function RMSFinfo = dataStreamRMSF(Ainfo)
%      --> Calculates RMSF on the DOF of the system based on the hinputted 
%          A-datamatrix that is in the format of the  dataMatrixInfo  data
%          structure. The output is formated according to the RMSFinfo
%          data structure. 
%
% function traitOut = traitProjection(fname,traitIn,U)
%      --> Projects statistical C-matrix into subspace defined by U basis 
%      --> vectors. Because the basis set changes, this is a rotation in 
%      --> the vector space, and, there is a dimension reduction to the 
%      --> size of the subspace spanned by the U basis vectors. 
%      --> Css = U' * C * U
%
% function traitOut = traitFilter(fname,traitIn,U,reduce)
%      --> Filters out from the C-matrices contained within traitIn the 
%      --> statistical correlations that are orthogonal to the U-subspace.
%      --> C_filtered = (U*U') * C * (U*U') 
%
% function Xshow = showDynamics(X)
% function Xshow = showDynamics(X,U)
%      --> From X trajectory the dynamics defined via the data matrix is
%      --> processed. Without specifying a basis set of vectors via U, the
%      --> function sparsifies the trajectory with an option to filter on 
%      --> U. If U defines a discriminant subspace, the motions that are
%      --> shown are functionally relevant. If U defines an indifference 
%      --> subspace, the motions show are conserved. 
%
% ========================================================================
%%                                                        write statements
version = 'version2020b';
end

