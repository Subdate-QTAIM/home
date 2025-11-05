#!/bin/bash

export KMP_STACKSIZE=2000M

path/to/Multiwfn "$1" << EOF | tee "${1%.wfn}.scr" #please, set the correct path to the Multiwfn binary
2
2
3
4
5
7
0
-4
4
6
0
-10
17
1
1
3
-4
3
7
2
1
-10
q
EOF
