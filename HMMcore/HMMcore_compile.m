% script with compilation commands

mex -setup % choose compiler options etc
mex('VBviterbi.c',['-' computer('arch')],'-O')
mex('VBviterbi_log.c',['-' computer('arch')],'-O')
