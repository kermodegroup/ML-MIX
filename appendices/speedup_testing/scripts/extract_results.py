import os
import re
import argparse

def find_loop_times(base_dir, output_file):
    # Prepare a list to store the results
    results = []

    # Walk through subdirectories and files
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".out"):
                file_path = os.path.join(root, file)
                subdirectory = os.path.relpath(root, base_dir)
                
                # Read the file and extract the second occurrence of "Loop time of"
                with open(file_path, 'r') as f:
                    content = f.read()
                    matches = re.findall(r"Loop time of (\d+\.\d+)", content)
                    
                    if len(matches) < 2:
                        print(f"error in {file_path}")
                    if len(matches) >= 2:
                        # Append the second match and the subdirectory
                        results.append((subdirectory, matches[1]))

    # Write the results to the output file
    with open(f"{args.base_directory}/{output_file}", 'w') as f:
        for subdirectory, loop_time in results:
            f.write(f"{subdirectory}\t{loop_time}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract loop times from .out files in a directory.")
    parser.add_argument("--base_directory","-bd", type=str, help="The base directory to search for .out files.")
    parser.add_argument("--output_file","-of", type=str, help="The output file to save results.")
    
    args = parser.parse_args()

    find_loop_times(args.base_directory, args.output_file)