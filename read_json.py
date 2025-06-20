import ijson
import json

def process(item):
    print(json.dumps({
        "features": {
            "source_num": int(item['features']['source_num']),
            "gene_num": int(item['features']['gene_num'])
        },
        "modality": str(item['modality']),
        "cell": str(item['cell']),
        "sample": int(item['sample']),
        "time_points": [int(tp) for tp in item['time_points']],
        "expression_rates": [float(rate) for rate in item['expression_rates']]
    }, indent=4))

filename = 'gene_expressions.json'

with open(filename, 'rb') as f:
    total = sum(1 for _ in ijson.items(f, 'item', use_float=True))

while True:
    print(f"\nValid range: 0 to {total}")
    try:
        start = int(input("Enter start index (inclusive): "))
        end = int(input("Enter end index (exclusive): "))
    except ValueError:
        print("Invalid input. Please enter integer values.")
        continue
    if not (0 <= start < end <= total):
        print(f"Invalid range. Please enter start >= 0, end <= {total}, and start < end.")
        continue

    with open(filename, 'rb') as f:
        for idx, item in enumerate(ijson.items(f, 'item', use_float=True)):
            if idx < start:
                continue
            if idx >= end:
                break
            process(item)

    another = input("\nWould you like to read another range? (y/n): ").lower()
    if another != 'y':
        break