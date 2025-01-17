#!/bin/bash
#SBATCH --job-name=download_urls
#SBATCH --output=download_%A_%a.log
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-6
#SBATCH --time=6:00:00

# Set the download directory
DOWNLOAD_DIR="/ocean/projects/cis240075p/asachan/datasets/SKM_multimodal_ageing"
mkdir -p "$DOWNLOAD_DIR"

# Set the input file containing URLs
URL_FILE="/ocean/projects/cis240075p/asachan/datasets/SKM_multimodal_ageing/bash_scripts/urls_cngb_processed_data.txt"

# Check if the URL file exists
if [ ! -f "$URL_FILE" ]; then
    echo "Error: URL file '$URL_FILE' not found!"
    exit 1
fi

# Create log files with array task ID
SUCCESS_LOG="successful_downloads_${SLURM_ARRAY_TASK_ID}.log"
FAILED_LOG="failed_downloads_${SLURM_ARRAY_TASK_ID}.log"
touch "$SUCCESS_LOG" "$FAILED_LOG"

# Create a temporary file with valid URLs only (no comments)
TEMP_URLS=$(mktemp)
# Only get lines that start with http, ignoring comments and empty lines
grep "^https" "$URL_FILE" > "$TEMP_URLS"

# Get total number of actual URLs (excluding comments and empty lines)
TOTAL_URLS=$(wc -l < "$TEMP_URLS")
echo "Total valid URLs found: $TOTAL_URLS"

# Verify we have URLs to process
if [ $TOTAL_URLS -eq 0 ]; then
    echo "Error: No valid URLs found in $URL_FILE"
    rm -f "$TEMP_URLS"
    exit 1
fi

# Calculate URLs per task (using actual number of URLs)
URLS_PER_TASK=3

# Calculate start and end indices for this task
START=$((SLURM_ARRAY_TASK_ID * URLS_PER_TASK + 1))
END=$((START + URLS_PER_TASK - 1))

echo "Task ${SLURM_ARRAY_TASK_ID} processing URLs $START to $END"
echo "Total URLs: $TOTAL_URLS"
echo "URLs per task: $URLS_PER_TASK"

# Process URLs for this task
sed -n "${START},${END}p" "$TEMP_URLS" | while read -r url; do
    # Skip if line is empty or doesn't start with https
    [[ -z "$url" || ! "$url" =~ ^https ]] && continue
    
    filename=$(basename "$url")
    echo "Downloading: $url"
    echo "Saving to: $DOWNLOAD_DIR/$filename"
    
    if wget --no-check-certificate \
         -O "$DOWNLOAD_DIR/$filename" \
         --tries=5 \
         --retry-connrefused \
         --waitretry=10 \
         --timeout=300 \
         --continue \
         "$url" 2>/dev/null; then
        echo "$(date): Successfully downloaded $url" >> "$SUCCESS_LOG"
        # Verify file size is non-zero
        if [ -s "$DOWNLOAD_DIR/$filename" ]; then
            echo "$(date): Verified file size for $filename" >> "$SUCCESS_LOG"
        else
            echo "$(date): Zero-size file detected for $url" >> "$FAILED_LOG"
            rm -f "$DOWNLOAD_DIR/$filename"  # Remove empty file
        fi
    else
        echo "$(date): Failed to download $url" >> "$FAILED_LOG"
    fi
    
    sleep 30
done

# Cleanup temporary file
rm -f "$TEMP_URLS"

echo "Task ${SLURM_ARRAY_TASK_ID} completed"

# If this is the last task, generate a summary
if [ "$SLURM_ARRAY_TASK_ID" -eq 6 ]; then
    echo "Generating final summary..."
    cat successful_downloads_*.log > final_successful_downloads.log
    cat failed_downloads_*.log > final_failed_downloads.log
    echo "Final Download Summary:" > download_summary.txt
    echo "Successfully downloaded: $(wc -l < final_successful_downloads.log)" >> download_summary.txt
    echo "Failed downloads: $(wc -l < final_failed_downloads.log)" >> download_summary.txt
fi