import json

def load(name):
    with open(f'{name}', "r") as f:
        data = json.load(f)
    
    return data