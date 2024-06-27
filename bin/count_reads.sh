#!/usr/bin/env bash
zcat -c $1 | wc -l | awk '{print $1/4}'
