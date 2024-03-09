#!/bin/bash

# Navigate to the NFDecoder directory (adjust the path as necessary)
cd NFDecoder

# Check if the AEResults directory exists, if not, create it
if [ ! -d "AEResults" ]; then
    echo "Creating AEResults directory..."
    mkdir AEResults
else
    echo "AEResults directory already exists."
fi


