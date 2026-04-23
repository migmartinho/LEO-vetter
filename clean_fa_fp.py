import csv

# Read FA/FP UIDs
fa_fp_uids = set()
with open('agg_fa_fp_tests_cleaned.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    for row in reader:
        fa_fp_uids.add(row[0])

# Read metrics and filter
metrics_rows = []
original_count = 0
with open('agg_metrics.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    metrics_rows.append(header)
    for row in reader:
        original_count += 1
        if row[0] in fa_fp_uids:
            metrics_rows.append(row)

# Write cleaned metrics
with open('agg_metrics_cleaned.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(metrics_rows)

print(f'Original metrics rows: {original_count}')
print(f'Cleaned metrics rows: {len(metrics_rows) - 1}')
print(f'Removed {original_count - (len(metrics_rows) - 1)} orphaned entries')