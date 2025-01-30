import os

def process_yaml_file(file_path):
    # Read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    updated_lines = []
    for line in lines:
        updated_lines.append(line)
        if "filename_core_out:" in line:
            # Capture the leading spaces
            leading_spaces = line[:len(line) - len(line.lstrip())]
            # Insert the new line with the same indentation
            updated_lines.append(f"{leading_spaces}datapath_out: Data/forecast\n")
    
    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)

# Process all YAML files in the directory
directory = "./"
for filename in os.listdir(directory):
    if filename.endswith('.yaml'):
        process_yaml_file(os.path.join(directory, filename))

print("Processing complete.")
