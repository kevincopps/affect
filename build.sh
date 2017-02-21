#!/usr/bin/env bash

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
    printf "\r\n"
}

filter_warnings() {
    # use awk to filter out the deprecated NumPy warning and 4 lines before it and one line after it
    #awk '/Using deprecated NumPy API/{for(x=NR-4;x<=NR+1;x++)d[x];}{a[NR]=$0}END{for(i=1;i<=NR;i++)if(!(i in d))print a[i]}' |
    perl -ne 'push @lines, $_;
          splice @lines, 0, 5 if /Using deprecated NumPy API/;
          print shift @lines if @lines > 4
          }{ print @lines;'
    #sed '/1 warning generated./d' |
    #sed '/cannot find cimported module/d'  # (caused with Python 3 cythonize)
}

echo "python setup.py build_ext -i"

python setup.py build_ext -i 2>&1 | filter_warnings && pipe_status

spinner


#while IFS= read -r result
#do
#    echo $result
#done < <(python setup.py build_ext -i 2>&1)
