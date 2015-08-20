#!/usr/bin/env python
import pandas as pd
import urllib2
from bs4 import BeautifulSoup
import re
import matplotlib.pyplot as plt
import cPickle as cpk
import sys

r = len(sys.argv)
