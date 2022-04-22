#!/bin/bash

MyProjectDir=`pwd`

echo $MyProjectDir
clang++ -I/usr/local/include -I$MyProjectDir/include SpheresScene.cpp -o CpuRayTracer