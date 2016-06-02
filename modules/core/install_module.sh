#!/bin/bash

INSTALL_DIR=$1

if [ -d "$INSTALL_DIR" ]; then

#I copy all the binaries
cp modules/core/bin/* $INSTALL_DIR/bin/

fi



