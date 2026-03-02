#!/bin/bash

source ./config

cd "$basedir/rep"

bwa index *.fa
