import csv
import sys

def filter_columns(input_file):
    # Open the input file (no output file argument)
    with open(input_file, mode='r', newline='', encoding='utf-8') as infile:
        # Set the delimiter to tab ('\t') for TSV files
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)  # Read the header (first row)
        
        # Identify columns that contain an underscore
        valid_columns = [col for col in header if '_' in col]
        
        # Now output the results to stdout (which will be redirected to a file)
        writer = csv.writer(sys.stdout, delimiter='\t')
        
        # Write the header of the filtered columns
        writer.writerow(valid_columns)
        
        # Loop through the rest of the rows and filter them
        for row in reader:
            # Create a new row with only the valid columns
            filtered_row = [row[i] for i in range(len(header)) if header[i] in valid_columns]
            writer.writerow(filtered_row)

if __name__ == "__main__":
    # Check if exactly one argument is provided for input file
    if len(sys.argv) != 2:
        print("Usage: python3 remove_columns.py input_file")
        sys.exit(1)

    # Get the input file path from the command-line argument
    input_file = sys.argv[1]

    # Call the function to filter the columns
    filter_columns(input_file)
