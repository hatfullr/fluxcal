#!/bin/bash

make clean

if [ $? -eq 0 ]; then
    make
    if [ $? -eq 0 ]; then
	echo ""
	echo "Success"
    else
	echo ""
	echo "Failed"
    fi
else
    echo ""
    echo "Failed"
fi
