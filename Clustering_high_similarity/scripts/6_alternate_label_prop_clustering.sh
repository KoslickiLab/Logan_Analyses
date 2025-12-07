#!/bin/bash
/usr/bin/time -v python manual_label_propogation.py --input-dir ../data/output/ --output-dir ../data/full_label_prop_viz_50k --verbose --max-plot-samples 50000
