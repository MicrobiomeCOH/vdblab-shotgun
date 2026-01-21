#!/usr/bin/env python3
import yaml
import glob
import os
import sys

def main(assembly_dir, config_file, output_file):
    assembly_dir = os.path.abspath(assembly_dir)
    # Load existing config
    with open(config_file) as f:
        config = yaml.safe_load(f)

    # Build assembly mapping
    assembly_mapping = {}
    for sample in config.get('R1', {}).keys():
        pattern = os.path.join(assembly_dir, f"megahit_{sample}.assembly.fasta")
        files = glob.glob(pattern)
        if files:
            assembly_mapping[sample] = files
        else:
            print(f"Warning: No assembly file found for sample {sample}")

    # Update config
    config['assembly'] = assembly_mapping

    # Save updated config to a new file
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    print(f"✅ Updated config written to {output_file}")
    print("Assembly section:")
    print(yaml.dump({'assembly': assembly_mapping}, default_flow_style=False))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python update_config.py <assembly_dir> <config_file> <output_file>")
        sys.exit(1)

    assembly_dir = sys.argv[1]
    config_file = sys.argv[2]
    output_file = sys.argv[3]
    main(assembly_dir, config_file, output_file)
