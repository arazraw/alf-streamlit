import csv
import re

def normalize_string(input_string):
    # Apply regex normalization
    input_string = re.sub(r'\b(?:of|the|and|och|Department|Dept|Dept.|Inst.|Institutionen|Inst)\b', '', input_string.lower())
    normalized_string = re.sub(r'[aeiouäöü\s.\-å]', '', input_string.lower())
    return normalized_string

input_file = 'affiliations_university.csv'
output_file = 'affiliations_university_norm.csv'

with open(input_file, 'r', newline='') as csvfile:
    with open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(csvfile)
        writer = csv.writer(outfile)
        for row in reader:
            if row:  # Check if row is not empty
                # Apply normalization to the first column
                normalized_string = normalize_string(row[0])
                # Write normalized string to the second column and copy the third column as is
                writer.writerow([row[0], normalized_string, row[2]])

print("Normalization complete. Output saved to:", output_file)
