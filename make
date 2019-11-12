#!/bin/bash

make clean

if [ $? -eq 0 ]; then
    make -j | tee build/compilation.output
    if [[ $? -ne 0 ]] && [[ "$(tail -n1 build/compilation.output)" != "<builtin>: recipe for target 'mathlib/odeint.o' failed" ]] && [[ "$(tail -n1 build/compilation.output)" != "<builtin>: recipe for target 'mathlib/rkqs.o' failed" ]]; then
	echo ""
	echo "Failed"
	exit 1
    fi
    
    while [[ "$(tail -n1 build/compilation.output)" == "<builtin>: recipe for target 'mathlib/odeint.o' failed" ]] || [[ "$(tail -n1 build/compilation.output)" == "<builtin>: recipe for target 'mathlib/rkqs.o' failed" ]]; do
	make -j | tee build/compilation.output
    done
    echo ""
    echo "Success"
    exit 0
else
    echo ""
    echo "Failed"
    exit 1
fi
