import json

def save_to_json(data: object, file_path: str) -> None:
    with open(f"{file_path}.json", "w") as json_file:
        json.dump(data, json_file)
    
def load_json(file_path: str) -> object:
    data = None
    with open(f"{file_path}.json", "r") as json_file:
        data = json.load(json_file)
    return data
