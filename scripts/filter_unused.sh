#!/bin/sh

# ==============================================================================
# Script Name: filter_unused.sh
# Description: Reads file paths from stdin and prints only those NOT currently
#              opened by any process (using lsof).
# Usage:       find . -name "*.txt" | ./filter_unused.sh
# ==============================================================================

# Check if lsof is installed
if ! command -v lsof &> /dev/null; then
    echo "Error: 'lsof' is not installed. Please install it to use this script." >&2
    exit 1
fi

# Read file paths line by line from standard input
# IFS= prevents trimming of leading/trailing whitespace
# -r prevents backslash interpretation
while IFS= read -r file; do
    
    # 1. Check if the file actually exists to avoid false positives/errors
    if [ -e "$file" ]; then
        
        # 2. Run lsof on the specific file.
        #    - lsof returns exit code 0 if the file IS open.
        #    - lsof returns exit code 1 if the file is NOT open.
        #    - We redirect stdout and stderr to /dev/null to keep output clean.
        if ! lsof "$file" >/dev/null 2>&1; then
            
            # 3. If lsof failed (file not in use), print the filename.
            echo "$file"
        fi
    fi
done
