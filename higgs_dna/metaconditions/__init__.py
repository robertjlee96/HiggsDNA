import os

metaconditions = {json_file.replace(".json", ""): json_file for json_file in [fl for fl in os.listdir(os.path.dirname(__file__)) if fl.endswith(".json")]}