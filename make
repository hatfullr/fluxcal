#!/bin/bash

make clean

if [ $? -eq 0 ]; then
    make
    if [ $? -eq 0 ]; then
	echo ""
	echo "Success"
	exit 0
    else
	echo ""
	echo "Failed"
	exit 1
    fi
else
    echo ""
    echo "Failed"
    exit 1
fi
