import os
import time
import argparse
import requests
import json
from Bio import SeqIO

# Base URL for the Dfam API
DFAM_BASE_URL = "https://dfam.org/api/"
# Endpoint for submitting new searches
SUBMIT_URL = f"{DFAM_BASE_URL}/searches"
# Template for checking search status and retrieving results
STATUS_URL_TEMPLATE = f"{DFAM_BASE_URL}/searches/{{search_id}}"

def submit_sequence(species, sequence, seq_id):
    """
    Submits a single sequence to the Dfam API for searching with fixed organism, cutoff, and evalue.

    Args:
        sequence (Bio.Seq.Seq): The sequence object from Biopython.
        seq_id (str): The identifier of the sequence.

    Returns:
        str: The search_id if submission is successful, None otherwise.
    """
    # Fixed values for organism, cutoff, and evalue as requeste
    fixed_organism = "Drosophila melanogaster"
    fixed_cutoff = "curated"
    fixed_evalue = "0.001" # Note: This will be sent even if cutoff is 'curated' as per the request

    # if starts with Bger
    if species.startswith("Bger") or species.startswith("Bgca") or species.startswith("Brac"):
        fixed_organism = "Halyomorpha halys"
    # if starts with Dsim
    elif species.startswith("D"):
        fixed_organism = "Drosophila melanogaster"
    # if starts with Phortica
    else:
        fixed_organism = "Other"

    payload = {
        "sequence": str(sequence),
        "name": seq_id,
        "organism": fixed_organism,
        "cutoff": fixed_cutoff,
        "evalue": fixed_evalue, # Always include this as per the request
    }

    try:
        response = requests.post(
            SUBMIT_URL,
            # Explicitly setting Content-Type to application/x-www-form-urlencoded
            headers={"Accept": "application/json", "Content-Type": "application/x-www-form-urlencoded"},
            data=payload, # requests.post with 'data' automatically handles form-urlencoded
            timeout=60
        )
        response.raise_for_status()
        # The API returns the search ID in the "id" field.
        # Successful response structure for initial submission:
        # {
        #   "id": "d60eb7f0-59ef-11f0-8c34-b9f7e1296e15",
        #   "code": 200
        # }
        return response.json()["id"]
    except requests.RequestException as e:
        print(f"❌ Error submitting {seq_id}: {e}")
        return None

def poll_search_result(search_id, max_attempts=30, wait_seconds=15):
    """
    Polls the Dfam API for the status and results of a submitted search.
    The polling response can now contain a 'status' field (e.g., 'RUNNING')
    and a 'message' field that might describe the running status or an error.
    Success is indicated by the presence of 'hits' or 'tandem_repeats' within the 'results' array.

    Args:
        search_id (str): The ID of the search to poll.
        max_attempts (int): Maximum number of attempts to poll.
        wait_seconds (int): Time to wait between attempts in seconds.

    Returns:
        dict: The search result dictionary if completed successfully, None if failed or timed out.
    """
    for attempt in range(max_attempts):
        try:
            response = requests.get(
                STATUS_URL_TEMPLATE.format(search_id=search_id),
                headers={"Accept": "application/json"},
                timeout=30
            )
            response.raise_for_status() # This will raise an exception for 4xx/5xx HTTP errors

            result = response.json()

            # Sample successful polling response (job completed):
            # {
            #   "submittedAt":"2025-07-05T22:34:34.000Z",
            #   "duration":"17 seconds",
            #   "searchParameters":"--species \"Drosophila melanogaster\" --cut_ga",
            #   "results":[
            #     {"query":"Input","length":"8833",
            #       "hits":[...],
            #       "tandem_repeats":[...]
            #     }
            #   ],
            #   "code":200
            # }

            # Sample polling response (job running/pending):
            # {
            #   "submittedAt":"2025-07-05T22:34:34.000Z",
            #   "duration":"Not finished",
            #   "searchParameters":"--species \"Drosophila melanogaster\" --cut_ga",
            #   "status":"RUNNING",
            #   "message":"Search running since Sat Jul 05 2025 15:34:42 GMT-0700 (Pacific Daylight Time).",
            #   "code":202
            # }

            # Sample polling response (job complete, preparing data - NOT an error):
            # {
            #   "submittedAt": "2025-07-05T22:30:06.000Z",
            #   "duration": "15 seconds",
            #   "searchParameters": "--species \"Drosophila melanogaster\" --cut_ga",
            #   "message": "Job Complete, Preparing Data",
            #   "code": 202
            # }

            # Sample error polling response:
            # {
            #   "submittedAt":"2025-07-05T22:36:00.000Z",
            #   "duration":"2 seconds",
            #   "searchParameters":"--species melanogaster --cut_ga",
            #   "status":"ERROR",
            #   "message":"Search failed."
            # }

            # First, check for successful completion (presence of 'results' array with 'hits' or 'tandem_repeats')
            if "results" in result and len(result["results"]) > 0 and \
               ("hits" in result["results"][0] or "tandem_repeats" in result["results"][0]):
                return result

            # Next, check for explicit status messages from the API
            if "status" in result:
                if result["status"] == "RUNNING" or result["status"] == "PEND":
                    # Search is still in progress
                    print(f"⏳ Waiting for search {search_id} to complete... Status: {result['status']} ({attempt+1}/{max_attempts})")
                    if "message" in result: # Include the message if it provides more context
                        print(f"   Message: {result['message']}")
                    time.sleep(wait_seconds)
                elif result["status"] == "ERROR" or result["status"] == "FAILED":
                    # Explicit error status
                    print(f"❌ Search {search_id} failed with status: {result['status']}")
                    if "message" in result:
                        print(f"   Details: {result['message']}")
                    return None
                else:
                    # Any other explicit status (not RUNNING/PENDING/ERROR/FAILED)
                    # For now, treat as still processing if not a known error status
                    print(f"⏳ Waiting for search {search_id} to complete... Status: {result.get('status', 'Unknown')} ({attempt+1}/{max_attempts})")
                    if "message" in result:
                        print(f"   Message: {result['message']}")
                    time.sleep(wait_seconds)
            elif "message" in result and result.get("code") == 202:
                # Specific case for "Job Complete, Preparing Data" with code 202, which is not an error
                print(f"⏳ Waiting for search {search_id} to complete... Message: {result['message']} (Code: {result['code']}) ({attempt+1}/{max_attempts})")
                time.sleep(wait_seconds)
                attempt-1
            elif "message" in result:
                # If there's a message and no specific status, and not the 202 "Job Complete" message,
                # it's likely an error message without a formal 'status' field.
                print(f"❌ Search {search_id} failed with message: {result['message']}")
                return None
            else:
                # If none of the above, still processing or unexpected response format
                print(f"⏳ Waiting for search {search_id} to complete... (No explicit status/results yet, attempt {attempt+1}/{max_attempts})")
                time.sleep(wait_seconds)

        except requests.RequestException as e:
            print(f"❌ Error polling {search_id}: {e}")
            time.sleep(wait_seconds) # Wait even on request error before retrying
    print(f"❌ Timed out waiting for search {search_id}")
    return None

def search_dfam(fasta_file, output_folder, delay=5):
    """
    Reads sequences from a FASTA file, submits them to Dfam, and saves results.

    Args:
        fasta_file (str): Path to the input FASTA file.
        output_folder (str): Path to the directory where results will be saved.
        delay (int): Delay in seconds between processing each sequence.
    """
    os.makedirs(output_folder, exist_ok=True)
    
    #use first 4 char of output folder to determine species/individual
    species = os.path.basename(output_folder)[:4]

    print(f"🔍 Reading sequences from {fasta_file}...")

    try:
        total_sequences = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at {fasta_file}")
        return
    except Exception as e:
        print(f"❌ Error parsing FASTA file {fasta_file}: {e}")
        return

    if total_sequences == 0:
        print("❌ No sequences found in FASTA. Exiting.")
        return

    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), 1):
        print(f"🔄 [{i}/{total_sequences}] Querying Dfam for: {record.id}")
        output_path = os.path.join(output_folder, f"{record.id}_dfam.json")

        if os.path.exists(output_path):
            print(f"⚠️ Output already exists for {record.id}, skipping...")
            continue

        if len(record.seq) == 0:
            warning = f"⚠️ Skipping {record.id}: Empty sequence.\n"
            with open(output_path, "w") as out_handle:
                out_handle.write(warning)
            print(warning.strip())
            continue

        if len(set(str(record.seq).upper())) == 1:
            warning = f"⚠️ Skipping {record.id}: Sequence consists of a single nucleotide.\n"
            with open(output_path, "w") as out_handle:
                out_handle.write(warning)
            print(warning.strip())
            continue

        # Call submit_sequence without organism, cutoff, evalue_string arguments
        search_id = submit_sequence(species, record.seq, record.id)
        if not search_id:
            continue

        result = poll_search_result(search_id)
        if result:
            with open(output_path, "w") as out_handle:
                json.dump(result, out_handle, indent=4)
            print(f"✅ Saved: {output_path}")
        else:
            print(f"❌ No result for {record.id}.")

        time.sleep(delay)

    print(f"✅ All Dfam queries completed for {fasta_file}.")


def main():
    """
    Main function to parse arguments and initiate the Dfam search.
    """
    parser = argparse.ArgumentParser(description="Query Dfam API with sequences from a FASTA file.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_folder", help="Output directory to save Dfam results")
    # Removed organism, cutoff, evalue arguments as they are now fixed in submit_sequence
    parser.add_argument("--delay", type=int, default=5, help="Delay between requests (default: 5 seconds)")

    args = parser.parse_args()

    search_dfam(
        fasta_file=args.fasta_file,
        output_folder=args.output_folder,
        delay=args.delay
    )


if __name__ == "__main__":
    main()
