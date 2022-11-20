import json


def load(path=None):
    if path is None:
        path = "../config.json"
    with open(path) as f:
        data = json.load(f)
    return data
