#!/bin/sh

acc importdata -p conf.json
acc calevent -p conf.json
acc migration -p conf.json
acc plotprofile -p conf.json
