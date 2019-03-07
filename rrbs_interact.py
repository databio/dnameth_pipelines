import argparse
import os
from pypiper import PipelineManager as PM

conf = os.path.expandvars(os.path.join("$CODE", "dnameth_pipelines", "src", "rrbs.yaml"))
assert os.path.isfile(conf), "Missing config file: {}".format(conf)
parser = argparse.ArgumentParser()
parser.add_argument("--config-file", required=True)
opts = parser.parse_args(["--config-file", conf])
pm = PM("RRBS", outfolder=os.getenv("HOME"), args=opts)
