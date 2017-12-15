#!/bin/bash

source ~/.bashrc
cd ~/bmxsim
HOSTNAME=`hostname -s`
celery -A bmxsim.celery_tasks worker -l info -c 10 &> celery_$HOSTNAME.log
