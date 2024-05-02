#!/bin/bash

MODEL_DESIGN="005-Mochi/model_design.txt"
OUTPUT_DIR="005-Mochi_output"
PROJECT_NAME="005-Mochi/mochi_model_pymochi_modabs"
HOLDOUT_ORDERS="2"
FEATURES="005-Mochi/features.txt"
CUSTOM_TRANSFORMATIONS="005-Mochi/custom_transformations.py"

run_mochi.py --model_design $MODEL_DESIGN --output_directory $OUTPUT_DIR --project_name $PROJECT_NAME --holdout_orders $HOLDOUT_ORDERS --features $FEATURES --custom_transformations $CUSTOM_TRANSFORMATIONS
005-Mochi/reformat_BindingMod.R