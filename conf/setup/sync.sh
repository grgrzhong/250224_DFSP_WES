#!/bin/bash

#############################################################################
## Sync local folders to a remote server using rsync
## Author: Zhong Guorui
## Date: 2025-04-15
#############################################################################

# Remote server details
remote_user="zhonggr"
remote_host="hpc" # hpc2021-io1.hku.hk 
remote_base_dir="/lustre1/g/path_my/250224_DFSP_WES"

# Define folders to sync (relative or absolute paths)
# Format: "local_path:remote_subpath"
folders_to_sync=(
    "data/reference:data/reference"
    "./containers:containers"
    # Add more folders as needed in the same format
)

# Function to ensure paths end with /
ensure_trailing_slash() {
    local path="$1"
    [[ "$path" != */ ]] && path="${path}/"
    echo "$path"
}

# Sync each folder
for folder_pair in "${folders_to_sync[@]}"; do
    # Split the pair into local and remote paths
    local_path="${folder_pair%%:*}"
    remote_subpath="${folder_pair#*:}"
    
    # Ensure paths have trailing slashes
    local_path=$(ensure_trailing_slash "$local_path")
    remote_subpath=$(ensure_trailing_slash "$remote_subpath")
    remote_base_dir=$(ensure_trailing_slash "$remote_base_dir")
    
    full_remote_path="${remote_base_dir}${remote_subpath}"
    
    echo "Syncing ${local_path} to ${remote_user}@${remote_host}:${full_remote_path}"
    
    # Create remote directory if it doesn't exist
    ssh ${remote_user}@${remote_host} "mkdir -p ${full_remote_path}"
    
    # Sync the folder
    rsync -av --update \
      "${local_path}" \
      "${remote_user}@${remote_host}:${full_remote_path}" \
      --dry-run
      
    echo "-----------------------------------"
done

# Remove --dry-run flag when you've confirmed it works correctly