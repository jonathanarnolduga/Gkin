# Gkin

To start the "solver", open a terminal session, cd to the GKIN directory and issue the following command:

java -cp server:classes theServer 8888 server/kin

To launch the GKIN GUI, open another terminal session, cd to the GKIN directory and issue the following command:

java -cp classes:server KinsolverGUI KinsolverGUI.properties

FYI, with the GKIN directory in /Applications as it is on the computer lab iMacs, here's the syntax for the two aliases I put in the users' login shell:

alias solverstart='java -cp /Applications/GKIN/server:/Applications/GKIN/classes theServer 8888 /Applications/GKIN/server/kin'

alias gkingui='java -cp /Applications/GKIN/classes:/Applications/GKIN/server KinsolverGUI /Applications/GKIN/KinsolverGUI.properties'

Of course, if the GKIN folder is somewhere other than /Applications, these lines will need to be edited accordingly.

DPB/JA
