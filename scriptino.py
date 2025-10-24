import json
import sys

# Check arguments
if len(sys.argv) != 5:
    print("Usage: python3 modify_n_tr.py <file_json> <nuovo_N> <nuovo_grad> <job_id>")
    sys.exit(1)

file_path = sys.argv[1]
new_n = int(sys.argv[2])
new_grad = float(sys.argv[3])
job_id = int(sys.argv[4])
new_tr = 1 + new_grad * new_n

# Function to modify `N`, `Tr` and `job_id`
def modify_n_tr_job(file_path, new_n, new_tr, job_id):
    try:
        # Read existing JSON
        with open(file_path, "r") as f:
            data = json.load(f)

        # Modify N
        if "iparams" in data and "N" in data["iparams"]:
            data["iparams"]["N"] = new_n
        else:
            print("Parameter 'N' not found in 'iparams'.")

        # Modify Tr
        if "dparams" in data and "Tr" in data["dparams"]:
            data["dparams"]["Tr"] = new_tr
        else:
            print("Parameter 'Tr' not found in 'dparams'.")

        # Modify job_id
        if "iparams" in data:
            data["iparams"]["job_id"] = job_id
        else:
            print("Parameter 'job_id' not found in 'iparams'.")

        # Save updated JSON
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)

        print(f"JSON updated successfully: N = {new_n}, Tr = {new_tr}, job_id = {job_id}")

    except Exception as e:
        print(f"Error modifying JSON file: {e}")
        sys.exit(1)

# Execute modification
modify_n_tr_job(file_path, new_n, new_tr, job_id)
