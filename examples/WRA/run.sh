#!/bin/sh

# import data from external
acc importdata -p conf.json

# calculate autocorrelograms
acc calevent -p conf.json

# migration for stations
acc migration -p conf.json

# plot profile
acc plotprofile -p conf.json
