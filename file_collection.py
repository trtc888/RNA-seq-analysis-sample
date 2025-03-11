## After downloading TSVs from TCGA, use this code to combine all files into 1 folder

import os
import shutil

def gather_tsv_files(source_dir, target_dir):
    
    # detect if the directories are exist
        if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Parse all subfolders
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith(".tsv"):
                
                # obtain full directories
                file_path = os.path.join(root, file)
                
                # Make full paths
                target_file_path = os.path.join(target_dir, file)
                
                # Change extension to avoid file replacement
                if os.path.exists(target_file_path):
                    base, extension = os.path.splitext(file)
                    counter = 1
                    while os.path.exists(target_file_path):
                        target_file_path = os.path.join(target_dir, f"{base}_{counter}{extension}")
                        counter += 1
                        
                # Copy files
                shutil.copy2(file_path, target_file_path)
                print(f"Copied {file_path} to {target_file_path}")

if __name__ == "__main__":
    
    # define directories
    source_directory = "/media/rt/code_data/Codes/rnaseq/TCGA_breast/gdc_data"
    target_directory = "/media/rt/code_data/Codes/rnaseq/TCGA_breast/TCGA_counts"

    # Run the code
    gather_tsv_files(source_directory, target_directory)
