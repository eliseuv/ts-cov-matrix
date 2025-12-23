import os
import re

path = "/run/media/evf/Research/contact_process/1d/all_active/error/"
params_dict = {}
for filename in os.listdir(path):
    # print(filename)
    pattern = r"(?:_|^)([a-zA-Z]\w*)=(\d+(?:\.\d+)?)"
    matches = re.findall(pattern, filename)
    params = dict(matches)
    params["diffusion"] = float(params["diffusion"])
    params["rate"] = float(params["rate"])
    print(params)
    # print((float(params["diffusion"]), float(params["rate"])))
