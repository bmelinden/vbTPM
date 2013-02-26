% s=VBviterbi_log(lnQ,lnqst) 
%
% Run the Viterbi algorithm using the log of the transition matrix (lnQ)
% and emission likelihood lnqst, e.g., find the most likely path for the
% following distribution 
% p(s(1:T)) ~  exp( sum_t  lnqst(t,s(t)) + lnQ(s(t-1),s(t)) ).
%
% This is intended for use with variation treatments of Hidden Markov
% Models, where the log of the path distribution comes out naturally from
% the compuation. The point of using th elog representation is to avoid
% numerical under- or overflow for models where the entries in lnQ and
% lnqst have very large magnitude.
%
% Martin Lindén, 2012-01-12

