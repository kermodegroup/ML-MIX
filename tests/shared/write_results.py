import os
import json
import re

def write_results(results, results_json, top_readme, insert_key, json_files):
    #check if test result json exists
    
    if os.path.exists(results_json):
        with open(results_json, "r") as f:
            data = json.load(f)
    else:
        data = {}

    #update data with new results
    data.update(results)

    #write data to json file
    with open(results_json, "w") as f:
        json.dump(data, f, indent=4)

    #now open both json files
    keys = ["Serial", "Parallel"]
    full_data = {}
    for i, json_file in enumerate(json_files):
        if os.path.exists(json_file):
            with open(json_file, "r") as f:
                full_data[keys[i]] = json.load(f)
        else:
            full_data[keys[i]] = {}

    #also write out markdown table in README.md
    with open("README.md", "w") as f:
        f.write("| Potential | Serial | Date and time | Parallel | Date and time |\n")
        f.write("| --- | --- | --- | --- | --- |\n")
        #for both test_types compile all keys
        all_keys = full_data["Serial"].keys() | full_data["Parallel"].keys()
        for key in all_keys:
            try:
                serial_result = full_data["Serial"][key]["result"]
                serial_datetime = full_data["Serial"][key]["datetime"]
            except KeyError:
                serial_result = "N/A"
                serial_datetime = "N/A"

            try:
                parallel_result = full_data["Parallel"][key]["result"]
                parallel_datetime = full_data["Parallel"][key]["datetime"]
            except KeyError:
                parallel_result = "N/A"
                parallel_datetime = "N/A"
            
            f.write(f"| {key} | {serial_result} | {serial_datetime} | {parallel_result} | {parallel_datetime} |\n")

    with open("README.md", "r") as f:
        table_content = f.read()

    #now open the top level read me
    with open(top_readme, "r") as f:
        readme_content = f.read()

    updated_readme_content = re.sub(
        fr"(<!-- {insert_key} start -->)(.*?)(<!-- {insert_key} end -->)",
        f"<!-- {insert_key} start -->\n{table_content}\n<!-- {insert_key} end -->",
        readme_content,
        flags=re.DOTALL,
    )

    with open(top_readme, "w") as f:
        f.write(updated_readme_content)