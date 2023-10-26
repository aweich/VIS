import os
import time
from pathlib import Path

configfile: "config.yml"
DATA = expand(config["data"]) 
PROCESS = os.path.join(config["processingdir"],str(config["experiment"]+"/")) #intermediate and results files are stored here

		
rule all:
	input:

rule Parsing_For_Vector:
	input:
	output:
	shell:
