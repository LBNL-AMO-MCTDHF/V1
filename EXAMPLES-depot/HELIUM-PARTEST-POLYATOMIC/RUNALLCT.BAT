#!/bin/bash

opt="Batchopt=2"

./RUN0.BAT.CT $opt
./RUN0.BAT.CTWRAP.DFFTPACK $opt
./RUN0.BAT.CT.DFFTPACK $opt
./RUN0.BAT.CTWRAP $opt

