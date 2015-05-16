#!/bin/bash

for opt in "" .vel .perp .perp.vel
do
./RunProjFlux.Bat 300 $opt
done

