#!/usr/bin/env bash

# specifying the CC= below may be need to avoid errors like "clang: warning: no such sysroot directory: '/Developer/SDKs/MacOSX10.6.sdk'"
##CC=/usr/bin/gcc
bold=`tput bold`
normal=`tput sgr0`

pipe_status() {
    if test ${PIPESTATUS[0]} -eq 0; then
        echo "${bold}build succeeded.${normal}"
    else
        echo "${bold}build failed.${normal}" 1>&2
        exit 1
    fi
}

spinner() {
    pid=$! # Process Id of the previous running command
    spin='-\|/'
    i=0
    while kill -0 $pid 2>/dev/null
    do
      i=$(( (i+1) %4 ))
      printf "\r${spin:$i:1}"
      sleep .1
    done
}

filter_numpy_deprecated() {
    # use awk to filter out the deprecated NumPy warning and 4 lines before it and one line after it
    awk '/Using deprecated NumPy API/{for(x=NR-4;x<=NR+1;x++)d[x];}{a[NR]=$0}END{for(i=1;i<=NR;i++)if(!(i in d))print a[i]}' |
    awk '!/1 warning generated./'
}

echo "python setup.py build_ext -i"
python setup.py build_ext -i 2>&1 | filter_numpy_deprecated && pipe_status &

spinner
