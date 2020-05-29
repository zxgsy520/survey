#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys

def printf(file):

    for line in open(file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        print(line)

printf(sys.argv[1])
