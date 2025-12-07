import os
from liftover import get_lifter

def convert_positions(input_path, output_path, chain_from, chain_to):
    print(f"Converting {input_path} from {chain_from} to {chain_to}...")
    
    # Initialize lifter (downloads chain file automatically if needed)
    converter = get_lifter(chain_from, chain_to)
    
    converted_count = 0
    skipped_count = 0
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                outfile.write(line) # Header or bad line
                continue
                
            chrom = parts[0] # usually 'chry' or 'chrY'
            # Ensure chrom format matches what liftover expects (usually 'chrY')
            if chrom.lower() == 'chry':
                liftover_chrom = 'chrY'
            else:
                liftover_chrom = chrom
                
            try:
                pos = int(parts[3])
                # LiftOver returns a list of matches, we take the first one
                new_coords = converter[liftover_chrom][pos]
                
                if new_coords and len(new_coords) > 0:
                    new_chrom, new_pos, strand = new_coords[0]
                    # Update position in the list
                    parts[3] = str(new_pos)
                    outfile.write('\t'.join(parts) + '\n')
                    converted_count += 1
                else:
                    # Could not map position (deleted in new assembly or complex region)
                    # print(f"Skipping {parts[1]}: {pos} not mapped")
                    skipped_count += 1
            except ValueError:
                outfile.write(line) # Probably header
            except Exception as e:
                print(f"Error on line: {line.strip()} -> {e}")
                skipped_count += 1

    print(f"Finished {input_path}. Converted: {converted_count}, Skipped/Unmapped: {skipped_count}")

# Paths
base_dir = "Yleaf/yleaf/data"
hg38_dir = os.path.join(base_dir, "hg38")
t2t_dir = os.path.join(base_dir, "t2t")

files_to_convert = [
    "new_positions.txt",
    "new_positions_ancient.txt", 
    "old_positions.txt",
    "old_positions_ancient.txt"
]

# Run conversion for all files
for filename in files_to_convert:
    in_file = os.path.join(hg38_dir, filename)
    out_file = os.path.join(t2t_dir, filename)
    
    if os.path.exists(in_file):
        # hg38 to hs1 (T2T-CHM13)
        convert_positions(in_file, out_file, 'hg38', 'hs1')
    else:
        print(f"Warning: {in_file} not found.")
