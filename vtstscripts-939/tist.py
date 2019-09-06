#!/usr/bin/env python
from __future__ import print_function 

import sys
import re
import os
import gzip

def sumIonicTime(f,subt,ast):
  "sum all ionic times in an OUTCAR"
  subtt = 0
  subett = 0
  elecStep = 0
  avgStepTime = 0 
  for line in f :
    string = line.split()
    if len(string) > 0 :
      if string[0] == 'LOOP+:' :
        temps = string[3]
        temps2 = temps[:-1]
        temps2 = float(temps2)
        subtt = subtt + temps2
      if string[0] == "LOOP:" :
        temps = string[3]
        temps2 = temps[:-1]
        temps2 = float(temps2)
        subett = subett + temps2
        elecStep = elecStep + 1
        avgStepTime = subett / elecStep
  return (subtt,avgStepTime)

def avgElectronicStep(f):
  "find the average electronic step time"
  subtt = 0 
  elecStep = 0
  avgStepTime = 0
  for line in f :
    print("aaa")
    string = line.split()
    print(string)
    if len(string) > 0 :
      if string[0] == "LOOP:" :
        temps = string[3]
        temps2 = temps[:-1]
        temps2 = float(temps2)
        subtt = subtt + temps2
        elecStep = elecStep + 1
        print(elecStep)
#  avgStepTime = subtt / elecStep
  return avgStepTime      
  
tt = 0
subtt = 0
avgElecStep = 0
aESSum = 0 # sum of avgElecStep
numOUTCAR = 0 # number of OUTCARS read
hasOUTCAR = os.path.isfile('OUTCAR')
if hasOUTCAR :
  numOUTCAR = numOUTCAR + 1
  print("There is a file OUTCAR :",hasOUTCAR)
  f = open("OUTCAR")
  (subtt,avgElecStep) = sumIonicTime(f,subtt,avgElecStep)
  print('this OUTCAR recorded total ionic step time:',subtt,' sec')
  tt = tt + subtt
  f.close()
  f = open("OUTCAR")
#  avgElecStep = avgElectronicStep(f)
  print("In this OUTCAR the average electronic step time is :", avgElecStep," sec")
  aESSum = aESSum + avgElecStep
#  f = open("OUTCAR")
else :
  print("There is no OUTCAR in current directory")
distList = os.listdir("./")
for dir in distList :
  if os.path.isdir(dir):
    print("current dir is ",dir)
    if os.path.isfile(os.path.join(dir,"OUTCAR.gz")) :
      numOUTCAR = numOUTCAR + 1
      print("reading from OUTCAR.gz")
      f = gzip.open(os.path.join(dir,"OUTCAR.gz"))
      (subtt,avgElecStep) = sumIonicTime(f,subtt,avgElecStep)
      print('this OUTCAR recorded total ionic step time:',subtt,' sec')
#      avgElecStep = avgElectronicStep(f)
      print("In this OUTCAR the average electronic step time is :", avgElecStep," sec")
      aESSum = aESSum + avgElecStep
      f.close()
      tt = tt + subtt
    else :
      print("There is no OUTCAR.gz here")
    if os.path.isfile(os.path.join(dir,"OUTCAR")) :
      numOUTCAR = numOUTCAR + 1
      subtt = 0
      f = open(os.path.join(dir,"OUTCAR"))
#      for line in f:
#        string = line.split()
#        if len(string) > 0 :
#          if string[0] == 'LOOP+:' :
#            temps = string[3]
#            temps2 = temps[:-1]
#            temps2 = float(temps2)
#            subtt = subtt + temps2
      (subtt,avgElecStep) = sumIonicTime(f,subtt,avgElecStep)
      print('this OUTCAR recorded total ionic step time:',subtt,' sec')
#      avgElecStep = avgElectronicStep(f)
      print("In this OUTCAR the average electronic step time is :", avgElecStep," sec")
      aESSum = aESSum + avgElecStep
      f.close()
      tt = tt + subtt
    else :
      print("There is no OUTCAR here")
    #print("yippi")
print("All OUTCARS combined took total ionic step time: ",tt,' sec')
print("which in hours is :",tt/3600,' hours')
print("checking for number of cores used...")
if os.path.isfile("job.sub") : 
  f = open("job.sub")
  for line in f:
    string = line.split()
    if len(string) > 1 :
      if string[1] == "-n" : 
        print("number of cores is", string[2])
        print("total computational time is :", tt/3600 * float(string[2])," hours")
  f.close() 
print("overall average electronic step time for all OUTCAR's read is :", aESSum / numOUTCAR," sec")


#for line in f:
   
