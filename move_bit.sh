#!/bin/bash
PROJECT_NAME="dy_cavity_relocker_2"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET_DIR="./pynq_overlays"
LATEST_DIR="$TARGET_DIR/latest"

mkdir -p "$TARGET_DIR"
mkdir -p "$LATEST_DIR"

echo "Monitoring Vivado bitstream generation for builds in $PROJECT_NAME..."
echo "Builds older than 14 days are automatically removed"
echo "Only bitstreams in the latest folder are kept for git pushes"

inotifywait -q -m -r --exclude "pynq_overlays" -e close_write --format '%w%f' ./build | while read NEW_BIT; do
	if [[ "$NEW_BIT" == *".bit" ]]; then 
		sleep 2
		TIMESTAMP=$(date +%Y%m%d_%H%M)
		BUILD_DIR="$TARGET_DIR/$TIMESTAMP"
		mkdir -p "$BUILD_DIR"


		LATEST_HWH=$(find . -path "./pynq_overlays" -prune -false -o -type f -name "*.hwh" -printf '%T+ %p\n' | sort -r | head -n 1 | cut -d' ' -f2-)


		if [ -f "$LATEST_HWH" ]; then
			cp "$NEW_BIT" "$BUILD_DIR/${PROJECT_NAME}.bit"
			cp "$LATEST_HWH" "$BUILD_DIR/${PROJECT_NAME}.hwh"


			cp "$NEW_BIT" "$LATEST_DIR/${PROJECT_NAME}.bit"
			cp "$LATEST_HWH" "$LATEST_DIR/${PROJECT_NAME}.hwh"



			echo "---------------------------------------------------------------"
			echo "PYNQ OVERLAY READY"
			echo "Archived locally in: $BUILD_DIR/"
			echo "Ready for Git in: $LATEST_DIR/"
			echo "---------------------------------------------------------------"

		else
			echo "WARNING: .bit saved, but no .hwh found"
		fi

		echo "Running auto-cleanup for builds older than 14 days"
		find "$TARGET_DIR" -mindepth 1 -maxdepth 1 -type d -name "20*" -mtime +14 -exec rm -rf {} +
	fi
done 
