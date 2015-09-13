#!/bin/bash

for opt in "" .vel .perp .perp.vel
do
./RunProjFlux.Bat 500 $opt
done

