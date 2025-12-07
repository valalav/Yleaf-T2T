import os
from liftover import get_lifter

def convert_positions(input_path, output_path, chain_from, chain_to):
    print(f"Converting {input_path} from {chain_from} to {chain_to}...")
    
    if not os.path.exists(input_path):
        print(f"Error: Input file not found: {input_path}")
        return

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
                
            chrom = parts[0] # usually 'chry'
            # LiftOver expects 'chrY' usually
            liftover_chrom = 'chrY'
                
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
                    # Could not map position
                    skipped_count += 1
            except ValueError:
                outfile.write(line) 
            except Exception as e:
                # print(f"Error on line: {line.strip()} -> {e}")
                skipped_count += 1

    print(f"Finished. Converted: {converted_count}, Skipped/Unmapped: {skipped_count}")

# Paths
base_dir = "Yleaf/yleaf/data"
hg38_dir = os.path.join(base_dir, "hg38")
t2t_dir = os.path.join(base_dir, "t2t")

# Ensure output dir exists
os.makedirs(t2t_dir, exist_ok=True)

# We only need to convert new_positions.txt as that's what the updater created
files_to_convert = ["new_positions.txt"]

for filename in files_to_convert:
    in_file = os.path.join(hg38_dir, filename)
    out_file = os.path.join(t2t_dir, filename)
    
    # hg38 to hs1 (T2T-CHM13)
    convert_positions(in_file, out_file, 'hg38', 'hs1')