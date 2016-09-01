#!/usr/bin/sh

cd jpas-z0p21
ls jpas-z0p21_*.txt > jpas-z0p21 && split jpas-z0p21 -nl/30 -d jpas-z0p21_
cd ..
cd jpas-z0p50
ls jpas-z0p50_*.txt > jpas-z0p50 && split jpas-z0p50 -nl/30 -d jpas-z0p50_
cd ..
cd  jpas-z0p90
ls jpas-z0p90_*.txt > jpas-z0p90 && split jpas-z0p90 -nl/30 -d jpas-z0p90_
cd ..
cd jplus-z0p21
ls jplus-z0p21_*.txt > jplus-z0p21 && split jplus-z0p21 -nl/30 -d jplus-z0p21_
cd ..
cd jplus-z0p50
ls jplus-z0p50_*.txt > jplus-z0p50 && split jplus-z0p50 -nl/30 -d jplus-z0p50_
cd ..
cd jplus-z0p90
ls jplus-z0p90_*.txt > jplus-z0p90 && split jplus-z0p90 -nl/30 -d jplus-z0p90_
cd ..

mkdir outs
mkdir outs/jpas-z0p21 outs/jpas-z0p50 outs/jpas-z0p90
mkdir outs/jplus-z0p21 outs/jplus-z0p50 outs/jplus-z0p90
mkdir outs/obs-jpas outs/obs-jplus
